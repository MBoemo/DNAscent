// -----------------------------------------------------------------------------
// Fused selective-scan custom op for the DNAscent Mamba model.
//
// Replaces the PyTorch `_chunked_ssm` path (Hillis-Steele parallel scan, which
// profiling showed is ~95% of GPU forward time and ~82% memory-traffic
// overhead) with a single fused, state-in-registers CUDA kernel.
//
// Registered as `dnascent::selective_scan` via TORCH_LIBRARY so the scripted
// TorchScript module (model.pt) can call it from both Python (export /
// validation) and the C++ LibTorch inference binary.
//
// Semantics (forward-time selective SSM, folding in the D-skip and z-gate):
//     a_t   = exp(dt[b,t,d] * A[d,:])                 (elementwise over n)
//     bx_t  = dt[b,t,d] * x[b,t,d] * B[b,t,:]         (elementwise over n)
//     h_t   = a_t (.) h_{t-1} + bx_t                  (h in R^N, per (b,d))
//     y_t   = ( sum_n h_t[n] * C[b,t,n] + x[b,t,d] * D[d] ) * silu(z[b,t,d])
//
// Bidirectionality and out_proj stay in Python; this op only does the
// forward-time scan for one direction.  Inference only: no autograd/backward.
// -----------------------------------------------------------------------------

#include <torch/library.h>
#include <ATen/ATen.h>
#include <c10/util/Exception.h>

#include <vector>
#include <cmath>

#if defined(DNASCENT_CUDA)
#include <ATen/cuda/CUDAContext.h>
#include <c10/cuda/CUDAGuard.h>
#include <cuda_runtime.h>

namespace {

constexpr int CHUNK   = 256;   // sequence chunk length (matches _chunked_ssm)
constexpr int MAX_N   = 32;    // upper bound on d_state held in registers (model uses 16)

__device__ __forceinline__ float siluf(float v) {
    return v / (1.0f + __expf(-v));
}

// Pass 1: per (b, chunk k, d) sequentially scan the chunk with zero initial
// state, emitting the composed affine map of the chunk:
//   carryA = product_t a_t        (coefficient on the incoming carry)
//   carryB = local h at chunk end (contribution assuming zero carry-in)
__global__ void ssm_pass1_kernel(
        const float* __restrict__ dt,   // (B,L,D)
        const float* __restrict__ A,    // (D,N)
        const float* __restrict__ Bm,   // (B,L,N)
        const float* __restrict__ x,    // (B,L,D)
        float* __restrict__ carryA,     // (B,NC,D,N)
        float* __restrict__ carryB,     // (B,NC,D,N)
        int B_, int L, int D, int N, int NC) {
    const int d = blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= D) return;
    const int k = blockIdx.y;
    const int b = blockIdx.z;

    const int t0 = k * CHUNK;
    if (t0 >= L) return;
    int t1 = t0 + CHUNK;
    if (t1 > L) t1 = L;

    float areg[MAX_N];
    float h[MAX_N];
    float aprod[MAX_N];
    for (int n = 0; n < N; ++n) {
        areg[n]  = A[d * N + n];
        h[n]     = 0.0f;
        aprod[n] = 1.0f;
    }

    for (int t = t0; t < t1; ++t) {
        const float dtv = dt[(b * L + t) * D + d];
        const float xv  = x[(b * L + t) * D + d];
        const int   boff = (b * L + t) * N;
        for (int n = 0; n < N; ++n) {
            const float a  = __expf(dtv * areg[n]);
            const float bx = dtv * xv * Bm[boff + n];
            h[n]     = a * h[n] + bx;
            aprod[n] = aprod[n] * a;
        }
    }

    const int cidx = ((b * NC + k) * D + d) * N;
    for (int n = 0; n < N; ++n) {
        carryA[cidx + n] = aprod[n];
        carryB[cidx + n] = h[n];
    }
}

// Carry scan: for each (b, d) sequentially compose the per-chunk affine maps to
// get the carry-in state at the start of every chunk.  NC is small (L/CHUNK),
// so this is cheap; parallelism is over B*D.
//   hin[0] = 0 ;  hin[k+1] = carryA[k] (.) hin[k] + carryB[k]
__global__ void ssm_carry_kernel(
        const float* __restrict__ carryA,  // (B,NC,D,N)
        const float* __restrict__ carryB,  // (B,NC,D,N)
        float* __restrict__ hin,           // (B,NC,D,N)
        int B_, int D, int N, int NC) {
    const int d = blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= D) return;
    const int b = blockIdx.y;

    float hc[MAX_N];
    for (int n = 0; n < N; ++n) hc[n] = 0.0f;

    for (int k = 0; k < NC; ++k) {
        const int idx = ((b * NC + k) * D + d) * N;
        for (int n = 0; n < N; ++n) {
            hin[idx + n] = hc[n];                                   // carry-in
            hc[n] = carryA[idx + n] * hc[n] + carryB[idx + n];      // advance
        }
    }
}

// Pass 2: re-scan each chunk starting from its carry-in state, recomputing the
// discretised coefficients (cheaper than storing the full (B,L,D,N) tensor),
// and write y with the D-skip and z-gate folded in.
__global__ void ssm_pass2_kernel(
        const float* __restrict__ dt,   // (B,L,D)
        const float* __restrict__ A,    // (D,N)
        const float* __restrict__ Bm,   // (B,L,N)
        const float* __restrict__ Cm,   // (B,L,N)
        const float* __restrict__ x,    // (B,L,D)
        const float* __restrict__ Dsk,  // (D,)
        const float* __restrict__ z,    // (B,L,D)
        const float* __restrict__ hin,  // (B,NC,D,N)
        float* __restrict__ y,          // (B,L,D)
        int B_, int L, int D, int N, int NC) {
    const int d = blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= D) return;
    const int k = blockIdx.y;
    const int b = blockIdx.z;

    const int t0 = k * CHUNK;
    if (t0 >= L) return;
    int t1 = t0 + CHUNK;
    if (t1 > L) t1 = L;

    float areg[MAX_N];
    float h[MAX_N];
    const int hidx = ((b * NC + k) * D + d) * N;
    for (int n = 0; n < N; ++n) {
        areg[n] = A[d * N + n];
        h[n]    = hin[hidx + n];
    }
    const float dsk = Dsk[d];

    for (int t = t0; t < t1; ++t) {
        const float dtv = dt[(b * L + t) * D + d];
        const float xv  = x[(b * L + t) * D + d];
        const int   boff = (b * L + t) * N;
        float acc = 0.0f;
        for (int n = 0; n < N; ++n) {
            const float a  = __expf(dtv * areg[n]);
            const float bx = dtv * xv * Bm[boff + n];
            h[n] = a * h[n] + bx;
            acc += h[n] * Cm[boff + n];
        }
        const float zv = z[(b * L + t) * D + d];
        y[(b * L + t) * D + d] = (acc + xv * dsk) * siluf(zv);
    }
}

}  // namespace

static at::Tensor selective_scan_cuda(
        const at::Tensor& dt, const at::Tensor& A, const at::Tensor& Bm,
        const at::Tensor& Cm, const at::Tensor& x, const at::Tensor& Dsk,
        const at::Tensor& z) {
    const at::Tensor dtc = dt.contiguous();
    const at::Tensor Ac  = A.contiguous();
    const at::Tensor Bc  = Bm.contiguous();
    const at::Tensor Cc  = Cm.contiguous();
    const at::Tensor xc  = x.contiguous();
    const at::Tensor Dc  = Dsk.contiguous();
    const at::Tensor zc  = z.contiguous();

    TORCH_CHECK(dtc.scalar_type() == at::kFloat, "selective_scan: dt must be float32");
    TORCH_CHECK(dtc.dim() == 3, "selective_scan: dt must be (B,L,D)");

    const int Bn = dtc.size(0);
    const int L  = dtc.size(1);
    const int D  = dtc.size(2);
    const int N  = Ac.size(1);
    TORCH_CHECK(N <= MAX_N, "selective_scan: d_state exceeds MAX_N");
    const int NC = (L + CHUNK - 1) / CHUNK;

    auto opts = dtc.options();
    at::Tensor carryA = at::empty({Bn, NC, D, N}, opts);
    at::Tensor carryB = at::empty({Bn, NC, D, N}, opts);
    at::Tensor hin    = at::empty({Bn, NC, D, N}, opts);
    at::Tensor y      = at::empty({Bn, L, D}, opts);

    const int threads = 128;
    const dim3 grid((D + threads - 1) / threads, NC, Bn);
    cudaStream_t stream = at::cuda::getCurrentCUDAStream();

    ssm_pass1_kernel<<<grid, threads, 0, stream>>>(
        dtc.data_ptr<float>(), Ac.data_ptr<float>(), Bc.data_ptr<float>(),
        xc.data_ptr<float>(), carryA.data_ptr<float>(), carryB.data_ptr<float>(),
        Bn, L, D, N, NC);

    const dim3 cgrid((D + threads - 1) / threads, Bn, 1);
    ssm_carry_kernel<<<cgrid, threads, 0, stream>>>(
        carryA.data_ptr<float>(), carryB.data_ptr<float>(), hin.data_ptr<float>(),
        Bn, D, N, NC);

    ssm_pass2_kernel<<<grid, threads, 0, stream>>>(
        dtc.data_ptr<float>(), Ac.data_ptr<float>(), Bc.data_ptr<float>(),
        Cc.data_ptr<float>(), xc.data_ptr<float>(), Dc.data_ptr<float>(),
        zc.data_ptr<float>(), hin.data_ptr<float>(), y.data_ptr<float>(),
        Bn, L, D, N, NC);

    return y;
}
#endif  // DNASCENT_CUDA

// -----------------------------------------------------------------------------
// CPU fallback: plain sequential recurrence.  Only used when no CUDA device is
// present (e.g. exporting/validating on a CPU box); correctness-first, not fast.
// -----------------------------------------------------------------------------
static at::Tensor selective_scan_cpu(
        const at::Tensor& dt, const at::Tensor& A, const at::Tensor& Bm,
        const at::Tensor& Cm, const at::Tensor& x, const at::Tensor& Dsk,
        const at::Tensor& z) {
    const at::Tensor dtc = dt.contiguous();
    const at::Tensor Ac  = A.contiguous();
    const at::Tensor Bc  = Bm.contiguous();
    const at::Tensor Cc  = Cm.contiguous();
    const at::Tensor xc  = x.contiguous();
    const at::Tensor Dc  = Dsk.contiguous();
    const at::Tensor zc  = z.contiguous();

    TORCH_CHECK(dtc.scalar_type() == at::kFloat, "selective_scan: dt must be float32");

    const int64_t Bn = dtc.size(0);
    const int64_t L  = dtc.size(1);
    const int64_t D  = dtc.size(2);
    const int64_t N  = Ac.size(1);

    at::Tensor y = at::empty({Bn, L, D}, dtc.options());

    const float* pdt = dtc.data_ptr<float>();
    const float* pA  = Ac.data_ptr<float>();
    const float* pB  = Bc.data_ptr<float>();
    const float* pC  = Cc.data_ptr<float>();
    const float* px  = xc.data_ptr<float>();
    const float* pD  = Dc.data_ptr<float>();
    const float* pz  = zc.data_ptr<float>();
    float* py = y.data_ptr<float>();

    std::vector<float> h(N);
    for (int64_t b = 0; b < Bn; ++b) {
        for (int64_t d = 0; d < D; ++d) {
            for (int64_t n = 0; n < N; ++n) h[n] = 0.0f;
            const float dsk = pD[d];
            for (int64_t t = 0; t < L; ++t) {
                const float dtv = pdt[(b * L + t) * D + d];
                const float xv  = px[(b * L + t) * D + d];
                const int64_t boff = (b * L + t) * N;
                float acc = 0.0f;
                for (int64_t n = 0; n < N; ++n) {
                    const float a  = std::exp(dtv * pA[d * N + n]);
                    const float bx = dtv * xv * pB[boff + n];
                    h[n] = a * h[n] + bx;
                    acc += h[n] * pC[boff + n];
                }
                const float zv = pz[(b * L + t) * D + d];
                const float silu = zv / (1.0f + std::exp(-zv));
                py[(b * L + t) * D + d] = (acc + xv * dsk) * silu;
            }
        }
    }
    return y;
}

// -----------------------------------------------------------------------------
// Registration.  The schema is visible to TorchScript scripting/loading; the
// backend-specific implementations are dispatched by device.
// -----------------------------------------------------------------------------
TORCH_LIBRARY(dnascent, m) {
    m.def("selective_scan(Tensor dt, Tensor A, Tensor B, Tensor C, Tensor x, "
          "Tensor D, Tensor z) -> Tensor");
}

TORCH_LIBRARY_IMPL(dnascent, CPU, m) {
    m.impl("selective_scan", selective_scan_cpu);
}

#if defined(DNASCENT_CUDA)
TORCH_LIBRARY_IMPL(dnascent, CUDA, m) {
    m.impl("selective_scan", selective_scan_cuda);
}
#endif
