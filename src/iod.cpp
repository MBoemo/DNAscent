//----------------------------------------------------------
// Copyright 2026 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

/*** Gillespie simulation of the following beacon calculus model:
 *
 * fast = 100000;
 * m = 2000;
 * v = 1.4;
 * fr = 0.00001;
 *
 * Fp[i] = {chr![i],fast}.[i <= m] -> {~chr?[i+1],v}.Fp[i+1];
 * Fm[i] = {chr![i],fast}.[i >= 0] -> {~chr?[i-1],v}.Fm[i-1];
 * Ori[i] = {fire,fr}.(Fm[i] || Fp[i]) + {chr?[i],fast};
 * SeedOris[i] = [i <= m] -> {seed, fast}.(Ori[i] || SeedOris[i+1]);
 *
 * SeedOris[0];
 *
 * Each kb position on the chromosome has an origin (Ori) that can
 * either fire at rate fr (spawning two replication forks) or be
 * passively replicated by an incoming fork (chr handshake).  Forks
 * move at rate v kb/min and terminate at chromosome boundaries or
 * when they encounter an already-replicated position.
 ***/

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <numeric>
#include <sstream>
#include <ctime>
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "htsInterface.h"
#include "iod.h"

#define CHR_LEN 1000 // Assumes reads of longer than 1 Mb are fleetingly rare and can be ignored for the purposes of IOD estimation
double MODEL_RES = 10.; // 100 bp resolution

static const char *help=
"meIODy: DNAscent executable that estimates inter-origin distance.\n"
"To run DNAscent meIODy, do:\n"
"   DNAscent meIODy -l /path/to/leftForks_DNAscent_forksense.bed \
                    -r /path/to/rightForks_DNAscent_forksense.bed \
					--origin /path/to/origins_DNAscent_forkSense.bed \
					--termination /path/to/terminations_DNAscent_forkSense.bed \
					-d /path/to/detectOutput.bam -o /path/to/output.IOD \
					--tPulse1 5. \
					--tPulse2 10. \
					-o /path/to/output.IOD \n"
"Required arguments are:\n"
"  -l,--left                 path to leftForks file from forkSense detect with `bed` extension,\n"
"  -r,--right                path to rightFork file from forkSense detect with `bed` extension,\n"
"     --origin               path to origin file from forkSense detect with `bed` extension,\n"
"     --termination          path to termination file from forkSense detect with `bed` extension,\n"
"  -d,--detect               path to output from detect with `detect` or `bam` extension,\n"
"     --tPulse1              duration (in minutes) of the first analogue pulse,\n"
"     --tPulse2              duration (in minutes) of the second analogue pulse,\n"
"  -o,--output               path to output file or directory for IOD.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


double wasserstein(std::vector<double> a, std::vector<double> b) {

	std::sort(a.begin(), a.end());
	std::sort(b.begin(), b.end());

	size_t n = a.size();
	size_t m = b.size();

	// 1-Wasserstein distance: integral of |F_a(x) - F_b(x)| dx
	double totalDist = 0.0;
	size_t i = 0, j = 0;
	double prevCDF_a = 0.0, prevCDF_b = 0.0;
	double prevVal = 0.0;
	bool first = true;

	while (i < n || j < m) {

		double val;
		if (i < n && (j >= m || a[i] <= b[j])) {
			val = a[i];
		} else {
			val = b[j];
		}

		if (!first) {
			double diff = prevCDF_a - prevCDF_b;
			totalDist += std::abs(diff) * (val - prevVal);
		}
		first = false;

		while (i < n && a[i] == val) {
			prevCDF_a += 1.0 / n;
			i++;
		}
		while (j < m && b[j] == val) {
			prevCDF_b += 1.0 / m;
			j++;
		}
		prevVal = val;
	}

	return totalDist;
}


double wassersteinPresorted(const std::vector<double> &a, const std::vector<double> &b) {

	size_t n = a.size();
	size_t m = b.size();

	double totalDist = 0.0;
	size_t i = 0, j = 0;
	double prevCDF_a = 0.0, prevCDF_b = 0.0;
	double prevVal = 0.0;
	bool first = true;

	while (i < n || j < m) {

		double val;
		if (i < n && (j >= m || a[i] <= b[j])) {
			val = a[i];
		} else {
			val = b[j];
		}

		if (!first) {
			double diff = prevCDF_a - prevCDF_b;
			totalDist += std::abs(diff) * (val - prevVal);
		}
		first = false;

		while (i < n && a[i] == val) {
			prevCDF_a += 1.0 / n;
			i++;
		}
		while (j < m && b[j] == val) {
			prevCDF_b += 1.0 / m;
			j++;
		}
		prevVal = val;
	}

	return totalDist;
}


struct Fork {
	int position;
	int direction; // +1 (plus strand) or -1 (minus strand)
};

struct ForkCall {
	int EdU_start;
	int EdU_end;
	int BrdU_start;
	int BrdU_end;
	int direction; // +1 (plus strand) or -1 (minus strand)
};

struct MatchedForkPair {
	ForkCall leftFork;
	ForkCall rightFork;
};


struct SimulationResult {
	std::vector<double> leftForkTimes;
	std::vector<double> rightForkTimes;
	std::vector<int> firedOriginPositions;
	std::vector<int> terminationPositions;
	double completionTime;
};


struct Arguments {
	std::string lForkInput;
	std::string rForkInput;
	std::string originInput;
	std::string terminationInput;
	std::string DetectInput;
	std::string output;
    bool specifiedLeft = false;
    bool specifiedRight = false;
    bool specifiedOrigin = false;
    bool specifiedTermination = false;
    bool humanReadable = false;
	bool specifiedOutput = false;
    bool specifiedDetect = false;
	bool specifiedPulse1 = false;
	bool specifiedPulse2 = false;
	int threads = 1;
	double EdUpulseLength; // in minutes
	double BrdUpulseLength; // in minutes
};


Arguments parseIODArguments( int argc, char** argv ) {

    if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent IOD." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}

	Arguments args;

	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-l" or flag == "--left" ){

 			if (i == argc-1) throw TrailingFlag(flag);

            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());

            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);

			args.lForkInput = strArg;
			i+=2;
			args.specifiedLeft = true;
		}
        else if ( flag == "-r" or flag == "--right" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.rForkInput = strArg;
			i+=2;
			args.specifiedRight = true;
        }
        else if ( flag == "--origin" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.originInput = strArg;
			i+=2;
			args.specifiedOrigin = true;
        }
		else if ( flag == "--termination" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.terminationInput = strArg;
			i+=2;
			args.specifiedTermination = true;
        }
		else if ( flag == "--tPulse1" ){

			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.EdUpulseLength = std::stod(strArg);
			i += 2; // Skip the flag and its argument
            args.specifiedPulse1 = true;
        }
        else if ( flag == "--tPulse2" ){
			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.BrdUpulseLength = std::stod(strArg);
			i += 2; // Skip the flag and its argument
            args.specifiedPulse2 = true;
        }
        else if ( flag == "-h" or flag == "--help" ){
            std::cout << help << std::endl;
            exit(EXIT_SUCCESS);
		}
        else if ( flag == "-d" or flag == "--detect" ){
 		
 			if (i == argc-1) throw TrailingFlag(flag);		
 		
 			std::string strArg( argv[ i + 1 ] );
 			
 			const char *ext = get_ext(strArg.c_str());
				
			if (strcmp(ext,"bam") == 0){
				args.humanReadable = false;
			}
			else if (strcmp(ext,"detect") == 0){
				args.humanReadable = true;
            }
			else{
				throw InvalidExtension(ext);
			}
            args.specifiedDetect = true;
			
			args.DetectInput = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){
		
			if (i == argc-1) throw TrailingFlag(flag);		
		
 			std::string strArg( argv[ i + 1 ] );
            args.output = strArg;
			i+=2;
			args.specifiedOutput = true;
        }
		else if (flag == "-t" or flag == "--threads") {

			if (i == argc-1) throw TrailingFlag(flag);
			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi(strArg);
			i += 2; // Skip the flag and its argument
		}
        else throw InvalidOption( flag );
	}
    if ( !args.specifiedOutput or !args.specifiedPulse1 or !args.specifiedPulse2 or !args.specifiedOrigin or !args.specifiedTermination or !args.specifiedDetect or (!args.specifiedLeft and !args.specifiedRight) ) {
        std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent IOD." << std::endl;
        exit(EXIT_FAILURE);
    }

    return args;
}


SimulationResult runGillespie(int m, double v, double fr, std::mt19937 &gen) {

	//map to model resolution
	m = m * MODEL_RES;
	v = v * MODEL_RES;
	fr = fr / MODEL_RES;

	// State: each position is either replicated or unreplicated (has a live Ori)
	std::vector<bool> replicated(m + 1, false);
	std::vector<Fork> forks;
	std::vector<int> firedOrigins;
	std::vector<int> terminationPositions;
	int nUnreplicated = m + 1;

	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	double t = 0.0;

    std::vector<double> leftForkTimes(m + 1, -1.0);
    std::vector<double> rightForkTimes(m + 1, -1.0);

	while (nUnreplicated > 0 || !forks.empty()) {

		int nForks = static_cast<int>(forks.size());

		// Propensities
		double fireRate = nUnreplicated * fr;  // any Ori can fire
		double forkRate = nForks * v;           // any fork can advance 1 kb
		double totalRate = fireRate + forkRate;

		if (totalRate <= 0.0) break;

		// Time to next event (exponentially distributed)
		double u1 = uniform(gen);
		if (u1 == 0.0) u1 = std::numeric_limits<double>::min();
		double dt = -std::log(u1) / totalRate;
		t += dt;

		// Select event
		double r = uniform(gen) * totalRate;

		if (r < fireRate) {

			// --- Origin firing event ---
			// Select the k-th unreplicated position
			int k = static_cast<int>(r / fr);
			if (k >= nUnreplicated) k = nUnreplicated - 1;

			int count = 0;
			int firePos = -1;
			for (int i = 0; i <= m; i++) {
				if (!replicated[i]) {
					if (count == k) {
						firePos = i;
						break;
					}
					count++;
				}
			}

			// Ori[firePos] fires: position becomes replicated, spawn Fm || Fp
			replicated[firePos] = true;
			nUnreplicated--;
			leftForkTimes[firePos] = t;
			rightForkTimes[firePos] = t;
			firedOrigins.push_back(firePos);

			forks.push_back({firePos, +1}); // Fp
			forks.push_back({firePos, -1}); // Fm

		} else {

			// --- Fork movement event ---
			int k = static_cast<int>((r - fireRate) / v);
			if (k >= nForks) k = nForks - 1;

			int newPos = forks[k].position + forks[k].direction;

			if (newPos < 0 || newPos > m || replicated[newPos]) {
				// Fork termination: boundary or collision
				forks[k] = forks.back();
				forks.pop_back();
				if (newPos >= 0 && newPos <= m && replicated[newPos]) {
					// Record termination position
					terminationPositions.push_back(newPos);
				}
			} else {
				// chr![newPos] syncs with Ori[newPos]: passive replication
				replicated[newPos] = true;
				nUnreplicated--;
				forks[k].position = newPos;
                if (forks[k].direction == -1) {
                    leftForkTimes[newPos] = t;
                } else {
                    rightForkTimes[newPos] = t;
                }
			}
		}
	}

	return {leftForkTimes, rightForkTimes, firedOrigins, terminationPositions, t};
}


std::pair<std::vector<ForkCall>, std::vector<ForkCall>> parseForkCalls(
	const SimulationResult &result,
	double pulseStart,
	double EdUpulseEnd,
	double pulseEnd,
	int readStart,
	int readEnd) {

	std::vector<ForkCall> leftForkCalls;
	std::vector<ForkCall> rightForkCalls;

	int n = std::min(readEnd + 1, static_cast<int>(result.leftForkTimes.size()));

	// Minimum length of EdU and BrdU tracks that would be realistically and reliably detected by the forkSense filters
	int BrdUReqLength = 2*MODEL_RES; 
	int EdUReqLength = 2*MODEL_RES;

	// Parse left fork tracks (direction = -1)
	// Contiguous segments where leftForkTimes >= 0 each correspond to one left fork
	int i = n - 1;
	while (readStart < i) {

		// Skip positions that weren't replicated by a left fork
		if (result.leftForkTimes[i] < 0) {
			i--;
			continue;
		}

		// Get contiguous segment of left fork replication (one fork track)
		int segStart = i;
		while (i >= readStart && result.leftForkTimes[i] >= 0) {
			i--;
		}
		int segEnd = i + 1;

		// Divide the segment into EdU and BrdU incorporation based on pulse lengths
		int eduStart = -1, eduEnd = -1, brduStart = -1, brduEnd = -1;
		for (int j = segStart; j >= segEnd; j--) {
			double t = result.leftForkTimes[j];
			if (t >= pulseStart && t < EdUpulseEnd) {
				if (eduStart == -1) eduStart = j;
				eduEnd = j;
			}
			if (t >= EdUpulseEnd && t < pulseEnd) {
				if (brduStart == -1) brduStart = j;
				brduEnd = j;
			}
		}
		

		// This length check ensures that this would actually be called as a fork by DNAscent
		int EdULength = (eduStart != -1 && eduEnd != -1) ? std::abs(eduEnd - eduStart + 1) : 0;
		int BrdULength = (brduStart != -1 && brduEnd != -1) ? std::abs(brduEnd - brduStart + 1) : 0;

		// check that BrdU track extends at least 2 kb from read start, otherwise this would likely be filtered out by DNAscent's read end QC
		if (eduStart != -1 && brduStart != -1 && EdULength >= EdUReqLength && BrdULength >= BrdUReqLength && (std::abs(eduStart - readEnd) >= 2*MODEL_RES || std::abs(brduEnd- readStart) >= 2*MODEL_RES)) { 
			ForkCall call;
			call.EdU_start = eduStart;
			call.EdU_end = eduEnd;
			call.BrdU_start = brduStart;
			call.BrdU_end = brduEnd;
			call.direction = -1;
			leftForkCalls.push_back(call);
		}
	}

	// Parse right fork tracks (direction = +1)
	n = std::min(readEnd + 1, static_cast<int>(result.rightForkTimes.size()));
	i = readStart;
	while (i < n) {

		// Skip positions that weren't replicated by a right fork
		if (result.rightForkTimes[i] < 0) {
			i++;
			continue;
		}

		// Get contiguous segment of right fork replication (one fork track)
		int segStart = i;
		while (i < n && result.rightForkTimes[i] >= 0) {
			i++;
		}
		int segEnd = i - 1;

		// Divide the segment into EdU and BrdU incorporation based on pulse lengths
		int eduStart = -1, eduEnd = -1, brduStart = -1, brduEnd = -1;
		for (int j = segStart; j <= segEnd; j++) {
			double t = result.rightForkTimes[j];
			if (t >= pulseStart && t < EdUpulseEnd) {
				if (eduStart == -1) eduStart = j;
				eduEnd = j;
			}
			if (t >= EdUpulseEnd && t < pulseEnd) {
				if (brduStart == -1) brduStart = j;
				brduEnd = j;
			}
		}
		
		// This length check ensures that this would actually be called as a fork by DNAscent
		int EdULength = (eduStart != -1 && eduEnd != -1) ? std::abs(eduEnd - eduStart + 1) : 0;
		int BrdULength = (brduStart != -1 && brduEnd != -1) ? std::abs(brduEnd - brduStart + 1) : 0;

		if (eduStart != -1 && brduStart != -1 && EdULength >= EdUReqLength && BrdULength >= BrdUReqLength && (std::abs(eduStart - readStart) >= 2*MODEL_RES || std::abs(brduEnd- readEnd) >= 2*MODEL_RES)) {
			ForkCall call;
			call.EdU_start = eduStart;
			call.EdU_end = eduEnd;
			call.BrdU_start = brduStart;
			call.BrdU_end = brduEnd;
			call.direction = +1;
			rightForkCalls.push_back(call);
		}
	}

	return {leftForkCalls, rightForkCalls};
}


std::pair<std::vector<double>, std::vector<double>> runPulseChase(SimulationResult &result, std::mt19937 &gen, std::vector<int> &readLengths, double EdUpulseTime, double BrdUpulseTime) {

	if (result.completionTime <= EdUpulseTime + BrdUpulseTime) return {};
	std::uniform_real_distribution<double> uniform(0.0, result.completionTime - EdUpulseTime - BrdUpulseTime);
	double pulseStart = uniform(gen);
	double EdUpulseEnd = pulseStart + EdUpulseTime;
	double pulseEnd = EdUpulseEnd + BrdUpulseTime;
	

    std::uniform_int_distribution<> readDist(0, readLengths.size()-1);
    int readIndex = readDist(gen);
	int midpoint = CHR_LEN / 2 * MODEL_RES;
	int sampledReadLength = readLengths[readIndex] * MODEL_RES; // in kb
	int readStart = midpoint - sampledReadLength / 2;
	int readEnd = midpoint + sampledReadLength / 2;
	if (readStart < 0 or readEnd >= CHR_LEN * MODEL_RES) {
		return {};
	}

	auto forkCalls = parseForkCalls(result, pulseStart, EdUpulseEnd, pulseEnd, readStart, readEnd);

	std::vector<ForkCall> &leftForks = forkCalls.first;
	std::vector<ForkCall> &rightForks = forkCalls.second;
	int totalForks = leftForks.size() + rightForks.size();

	if (totalForks != 1) return {};

	const ForkCall &fork = leftForks.empty() ? rightForks[0] : leftForks[0];
	double behindDist, aheadDist;
	if (fork.direction == +1) {
		behindDist = static_cast<double>( (fork.EdU_start - readStart) / MODEL_RES );
		aheadDist = static_cast<double>( (readEnd - fork.BrdU_end) / MODEL_RES );
	} else {
		behindDist = static_cast<double>( (readEnd - fork.EdU_start) / MODEL_RES );
		aheadDist = static_cast<double>( (fork.BrdU_end - readStart) / MODEL_RES );
	}

	return {{behindDist}, {aheadDist}};
}


std::vector<double> parseForkBed(std::string fileInput, std::vector<std::string> &ignoreIDs) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
	std::vector<double> forkLengths;
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::vector< std::string > columns = split(line);               
            std::string readID = columns[3];

            if ( std::find(ignoreIDs.begin(), ignoreIDs.end(), readID) == ignoreIDs.end() ) {

                int pulse5Prime = std::stoi(columns[1]);
                int pulse3Prime = std::stoi(columns[2]);
                int read5Prime = std::stoi(columns[4]);
                int read3Prime = std::stoi(columns[5]);

				assert(read5Prime <= pulse3Prime && pulse3Prime <= read3Prime);
				assert(read5Prime <= pulse5Prime && pulse5Prime <= read3Prime);

				// ignore forks near the end of the read
				if (std::abs(pulse5Prime - read5Prime) < 3000) continue; 
				if (std::abs(read3Prime - pulse3Prime) < 3000) continue; 
				if (std::abs(pulse5Prime - read3Prime) < 3000) continue; 
				if (std::abs(pulse3Prime - read5Prime) < 3000) continue; 

				forkLengths.push_back( (pulse3Prime - pulse5Prime)/1000.); // in kb
			}
		}
	}

    file.close();
    return forkLengths;
}


void getIgnoreIDs(std::string fileInput, std::vector<std::string> &ignoreIDs) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::vector< std::string > columns = split(line);               
            std::string readID = columns[3];

			 if ( std::find(ignoreIDs.begin(), ignoreIDs.end(), readID) == ignoreIDs.end() ) {

				 ignoreIDs.push_back(readID);
			}
		}
	}
    file.close();
}


std::pair<std::vector<double>, std::vector<double>> calcSingleForkDistances(
	const std::string &leftForkFile,
	const std::string &rightForkFile) {

	struct ForkInfo {
		int pulse5Prime;
		int pulse3Prime;
		int direction;  // +1 (rightward) or -1 (leftward)
		int readStart;  // 5' read boundary
		int readEnd;    // 3' read boundary
	};

	// Group forks by readID
	std::map<std::string, std::vector<ForkInfo>> readForks;

	// Parse left fork BED
	{
		std::ifstream file(leftForkFile);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << leftForkFile << "\n";
			exit(EXIT_FAILURE);
		}
		std::string line;
		while (std::getline(file, line)) {
			if (!line.empty() && line[0] != '#') {
				std::vector<std::string> columns = split(line);
				std::string readID = columns[3];
				int pulse5Prime = std::stoi(columns[1]);
				int pulse3Prime = std::stoi(columns[2]);
				int read5Prime = std::stoi(columns[4]);
				int read3Prime = std::stoi(columns[5]);

				// We're taking the max distance anyway, so both of these need to fail to throw a read out
				if (std::abs(pulse3Prime - read3Prime) < 2000 and std::abs(pulse5Prime - read5Prime) < 2000) continue;

				readForks[readID].push_back({pulse5Prime, pulse3Prime, -1, read5Prime, read3Prime});
			}
		}
	}

	// Parse right fork BED
	{
		std::ifstream file(rightForkFile);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << rightForkFile << "\n";
			exit(EXIT_FAILURE);
		}
		std::string line;
		while (std::getline(file, line)) {
			if (!line.empty() && line[0] != '#') {
				std::vector<std::string> columns = split(line);
				std::string readID = columns[3];
				int pulse5Prime = std::stoi(columns[1]);
				int pulse3Prime = std::stoi(columns[2]);
				int read5Prime = std::stoi(columns[4]);
				int read3Prime = std::stoi(columns[5]);

				// We're taking the max distance anyway, so both of these need to fail to throw a read out
				if (std::abs(pulse3Prime - read3Prime) < 2000 and std::abs(pulse5Prime - read5Prime) < 2000) continue;

				readForks[readID].push_back({pulse5Prime, pulse3Prime, +1, read5Prime, read3Prime});
			}
		}
	}

	// Only process reads with exactly one fork call
	std::vector<double> behindDists;
	std::vector<double> aheadDists;
	for (const auto &entry : readForks) {
		const std::vector<ForkInfo> &forks = entry.second;
		if (forks.size() != 1) continue;

		const ForkInfo &fork = forks[0];
		if (fork.direction == +1) {
			// Rightward fork: behind = EdU start to 5' end, ahead = BrdU end to 3' end
			behindDists.push_back(static_cast<double>(fork.pulse5Prime - fork.readStart) / 1000.0);
			aheadDists.push_back(static_cast<double>(fork.readEnd - fork.pulse3Prime) / 1000.0);
		} else {
			// Leftward fork: behind = EdU start to 3' end, ahead = BrdU end to 5' end
			behindDists.push_back(static_cast<double>(fork.readEnd - fork.pulse3Prime) / 1000.0);
			aheadDists.push_back(static_cast<double>(fork.pulse5Prime - fork.readStart) / 1000.0);
		}
	}

	return {behindDists, aheadDists};
}


void detectReadlength(std::string filename, std::vector<int> &readLengths) {

    std::ifstream detectFile(filename);
    std::string line;
    int prog = 0;

	while( std::getline( detectFile, line ) ){

        if (line.substr(0,1) == "#" or line.length() == 0) continue; //ignore header and blank lines
        if ( line.substr(0,1) == ">" ){

            std::vector< std::string > columns = split(line);
            assert(columns.size() == 5); //well-formed detect header

            int refStart = std::stoi(columns[2]);
            int refEnd = std::stoi(columns[3]);
            
            readLengths.push_back((refEnd - refStart)/1000); // convert to kb

            prog ++;
            if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
        }
    }
	std::cout << " Done." << std::endl;
    detectFile.close();
}


void bamReadlength(std::string filename, std::vector<int> &readLengths) {

    htsFile *bam_fh = sam_open(filename.c_str(), "r");
    if (bam_fh == NULL) throw IOerror(filename);

    bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
    bam1_t *itr_record = bam_init1();

    int prog = 0;

    while(sam_read1(bam_fh, bam_hdr, itr_record) >= 0){
          
        int refStart,refEnd;
        getRefEnd(itr_record,refStart,refEnd);
        readLengths.push_back((refEnd - refStart)/1000); // convert to kb

        prog ++;
        if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
    }
	std::cout << " Done." << std::endl;
	bam_destroy1(itr_record);				
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
}


int iod_main(int argc, char** argv) {

	Arguments args = parseIODArguments(argc, argv);

	// Get readIDs to ignore (those with an origin or termination call) so that we get an accurate measure of fork speed
	std::vector<std::string> ignoreIDs;
	getIgnoreIDs(args.originInput, ignoreIDs);
	getIgnoreIDs(args.terminationInput, ignoreIDs);

	// Calculate median fork speed from bed files, ignoring any forks near ends of reads or reads with an origin or termination call
	std::vector<double> forkLengths = parseForkBed(args.lForkInput, ignoreIDs);
	std::vector<double> rightForkLengths = parseForkBed(args.rForkInput, ignoreIDs);
	forkLengths.insert(forkLengths.end(), rightForkLengths.begin(), rightForkLengths.end());
	double avgForkSpeed = vectorMean(forkLengths) / (args.EdUpulseLength + args.BrdUpulseLength); // convert to kb/min
	std::cout << "Mean fork speed = " << std::fixed << std::setprecision(3) << avgForkSpeed << " kb/min" << std::endl;

	// Calculate behind and ahead distances for single-fork reads
	auto dataDists = calcSingleForkDistances(args.lForkInput, args.rForkInput);
	std::vector<double> data_behindDists = dataDists.first;
	std::vector<double> data_aheadDists = dataDists.second;
	std::cout << "Total distance measurements: " << data_behindDists.size() << std::endl;

	// Compute max(behind, ahead) for each read — the less-truncated flank
	std::vector<double> data_maxDists(data_behindDists.size());
	for (size_t i = 0; i < data_behindDists.size(); i++) {
		data_maxDists[i] = std::max(data_behindDists[i], data_aheadDists[i]);
	}

	// Precompute mean for normalising the Wasserstein component
	double meanMax = vectorMean(data_maxDists);

    // Parse detect output
    std::vector<int> readLengths;

    if (args.humanReadable) {// is .detect format

        detectReadlength(args.DetectInput, readLengths);
    }
    else{// is modbam
        
        bamReadlength(args.DetectInput, readLengths);       
    }

	int nSims = 20000; 
	unsigned int seed = 42;
	std::pair<std::vector<double>, std::vector<double>> *saveSimDist = nullptr;
	std::vector<std::tuple<double, double, double>> landscapePoints; // (fr, IOD, W)
	const size_t MAX_LANDSCAPE_SAMPLES = 50000;

	struct LandscapeEntry {
		std::vector<double> sortedMax;
		std::vector<double> sortedBehind;
		std::vector<double> sortedAhead;
		double medianIOD;
	};
	std::map<double, LandscapeEntry> landscapeSimDists; // fr -> stored sim distributions + median IOD

	// Objective function: run simulations at a given fr and return combined Wasserstein distance to data
	auto objective = [&](double fr, double &medianIOD) -> double {
		std::vector<double> sim_behindDists;
		std::vector<double> sim_aheadDists;
		std::vector<int> sim_IODs;

		#pragma omp parallel for num_threads(args.threads) schedule(dynamic)
		for (int sim = 0; sim < nSims; sim++) {

			std::mt19937 gen(seed + sim);
			SimulationResult result = runGillespie(CHR_LEN, avgForkSpeed, fr, gen);
			std::vector<double> behindDists;
			std::vector<double> aheadDists;
			std::vector<int> iods;
			
			for (size_t i = 0; i < 100; i++){
			
				auto distances = runPulseChase(result, gen, readLengths, args.EdUpulseLength, args.BrdUpulseLength);
				behindDists.insert(behindDists.end(), distances.first.begin(), distances.first.end());
				aheadDists.insert(aheadDists.end(), distances.second.begin(), distances.second.end());

				// Compute IODs from fired origin positions
				std::vector<int> origins = result.firedOriginPositions;
				std::sort(origins.begin(), origins.end());
				for (size_t j = 1; j < origins.size(); j++) {
					iods.push_back( (origins[j] - origins[j - 1]) / MODEL_RES ); // convert to kb
				}
			}

			#pragma omp critical
			{
				sim_behindDists.insert(sim_behindDists.end(), behindDists.begin(), behindDists.end());
				sim_aheadDists.insert(sim_aheadDists.end(), aheadDists.begin(), aheadDists.end());
				sim_IODs.insert(sim_IODs.end(), iods.begin(), iods.end());
			}
		}
		medianIOD = vectorMedian(sim_IODs);

		// Compute max(behind, ahead) for each simulated read
		std::vector<double> sim_maxDists(sim_behindDists.size());
		for (size_t i = 0; i < sim_behindDists.size(); i++) {
			sim_maxDists[i] = std::max(sim_behindDists[i], sim_aheadDists[i]);
		}

		// Store subsampled, sorted sim distributions for bootstrap CI
		{
			auto subsample = [](const std::vector<double> &values, size_t maxN) -> std::vector<double> {
				std::vector<double> out;
				if (values.size() > maxN) {
					out.resize(maxN);
					std::mt19937 subGen(42);
					std::uniform_int_distribution<size_t> idx(0, values.size() - 1);
					for (size_t s = 0; s < maxN; s++) {
						out[s] = values[idx(subGen)];
					}
				} else {
					out = values;
				}
				std::sort(out.begin(), out.end());
				return out;
			};

			LandscapeEntry entry;
			entry.sortedMax = subsample(sim_maxDists, MAX_LANDSCAPE_SAMPLES);
			entry.sortedBehind = subsample(sim_behindDists, MAX_LANDSCAPE_SAMPLES);
			entry.sortedAhead = subsample(sim_aheadDists, MAX_LANDSCAPE_SAMPLES);
			entry.medianIOD = medianIOD;
			landscapeSimDists[fr] = std::move(entry);
		}

		if (saveSimDist) *saveSimDist = {sim_behindDists, sim_aheadDists};

		// Objective: Wasserstein distance on max(behind, ahead)
		double w = (!sim_maxDists.empty() && !data_maxDists.empty()) ? wasserstein(sim_maxDists, data_maxDists) / meanMax : 0.0;
		landscapePoints.push_back(std::make_tuple(fr, medianIOD, w));
		return w;
	};

	// Grid search in log10(fr) space: coarse pass across the full range,
	// then a fine pass around the minimum.
	double logFrMin = std::log10(1e-5);
	double logFrMax = std::log10(1e-2);

	std::cerr << "Coarse grid search for optimal firing rate..." << std::endl;
	const int nCoarse = 20;
	double coarseBestW = std::numeric_limits<double>::max();
	double coarseBestLogFr = logFrMin;
	for (int i = 0; i < nCoarse; i++) {
		double logFr = logFrMin + i * (logFrMax - logFrMin) / (nCoarse - 1);
		double fr = std::pow(10.0, logFr);
		double iod;
		double w = objective(fr, iod);
		std::cerr << "[Coarse " << std::setw(2) << (i + 1) << "/" << nCoarse << "] "
				  << "fr = " << std::scientific << std::setprecision(3) << std::setw(9) << fr
				  << "  W = " << std::fixed << std::setprecision(3) << std::setw(5) << w << std::endl;
		if (w < coarseBestW) {
			coarseBestW = w;
			coarseBestLogFr = logFr;
		}
	}

	std::cerr << "Fine grid search around coarse minimum..." << std::endl;
	const int nFine = 20;
	double fineLogLow = std::max(logFrMin, coarseBestLogFr - 0.5);
	double fineLogHigh = std::min(logFrMax, coarseBestLogFr + 0.5);
	for (int i = 0; i < nFine; i++) {
		double logFr = fineLogLow + i * (fineLogHigh - fineLogLow) / (nFine - 1);
		// Skip if we already have a nearby point
		bool exists = false;
		for (const auto &lp : landscapeSimDists) {
			if (std::abs(std::log10(lp.first) - logFr) < 0.02) { exists = true; break; }
		}
		if (exists) continue;
		double fr = std::pow(10.0, logFr);
		double iod;
		double w = objective(fr, iod);
		std::cerr << "[Fine  " << std::setw(2) << (i + 1) << "/" << nFine << "] "
				  << "fr = " << std::scientific << std::setprecision(3) << std::setw(9) << fr
				  << "  W = " << std::fixed << std::setprecision(3) << std::setw(5) << w << std::endl;
	}

	// Determine the point estimate from the landscape: sort data once and
	// scan all landscape entries with the same procedure the bootstrap uses,
	// so the point estimate and CI are guaranteed to be consistent.
	std::vector<double> data_max_sorted = data_maxDists;
	std::sort(data_max_sorted.begin(), data_max_sorted.end());

	double bestFr = 0.0;
	double bestW = std::numeric_limits<double>::max();
	double bestIOD = 0.0;
	for (const auto &lp : landscapeSimDists) {
		double w = wassersteinPresorted(lp.second.sortedMax, data_max_sorted) / meanMax;
		if (w < bestW) {
			bestW = w;
			bestFr = lp.first;
			bestIOD = lp.second.medianIOD;
		}
	}
	std::cerr << "Optimal: fr = " << std::scientific << std::setprecision(3) << bestFr
			  << "  IOD = " << std::fixed << std::setprecision(0) << bestIOD
			  << "  W = " << std::fixed << std::setprecision(3) << bestW << std::endl;

	// Evaluate at the optimum to get sim distributions for output
	std::pair<std::vector<double>, std::vector<double>> bestSimDist;
	saveSimDist = &bestSimDist;
	double iodFinal;
	objective(bestFr, iodFinal);
	saveSimDist = nullptr;

	// Bootstrap CI: resample the observed data and find the best landscape point for each resample
	const int nBootstrap = 10000;
	std::vector<double> bootstrapIODs;
	bootstrapIODs.reserve(nBootstrap);

	std::mt19937 bootGen(42);

	// Helper to compute smoothed bootstrap parameters for a distribution
	struct SmoothParams {
		double mean;
		double smoothBW;
		double shrinkFactor;
	};
	auto computeSmoothParams = [](const std::vector<double> &data) -> SmoothParams {
		SmoothParams p;
		p.mean = 0.0;
		for (const auto &v : data) p.mean += v;
		p.mean /= data.size();
		double sd = 0.0;
		for (const auto &v : data) sd += (v - p.mean) * (v - p.mean);
		sd = std::sqrt(sd / (data.size() - 1));
		p.smoothBW = sd * std::pow(static_cast<double>(data.size()), -0.2) * 0.5;
		p.shrinkFactor = 1.0 / std::sqrt(1.0 + (p.smoothBW * p.smoothBW) / (sd * sd));
		return p;
	};

	SmoothParams spMax = computeSmoothParams(data_maxDists);

	for (int b = 0; b < nBootstrap; b++) {

		// Resample (behind, ahead) pairs together with replacement and smooth
		size_t nData = data_behindDists.size();
		std::uniform_int_distribution<size_t> readIdx(0, nData - 1);
		std::normal_distribution<double> noiseMax(0.0, spMax.smoothBW);

		std::vector<double> bootMax(nData);
		for (size_t i = 0; i < nData; i++) {
			size_t idx = readIdx(bootGen);
			double mval = data_maxDists[idx] + noiseMax(bootGen);
			mval = std::max(0.0, spMax.mean + (mval - spMax.mean) * spMax.shrinkFactor);
			bootMax[i] = mval;
		}
		std::sort(bootMax.begin(), bootMax.end());

		// Compute combined Wasserstein distance at each landscape point for this bootstrap sample
		std::vector<std::pair<double, double>> bootLandscape; // (log10(fr), W)
		std::vector<std::pair<double, double>> iodLandscape;  // (log10(fr), IOD)
		for (const auto& lp : landscapeSimDists) {
			double logFr = std::log10(lp.first);
			double w = wassersteinPresorted(lp.second.sortedMax, bootMax) / meanMax;
			bootLandscape.push_back({logFr, w});
			iodLandscape.push_back({logFr, lp.second.medianIOD});
		}

		// Find the index of the minimum W
		size_t minIdx = 0;
		for (size_t i = 1; i < bootLandscape.size(); i++) {
			if (bootLandscape[i].second < bootLandscape[minIdx].second) minIdx = i;
		}

		// Fit a quadratic to the minimum and its neighbours to interpolate
		double bootIOD;
		if (minIdx > 0 && minIdx < bootLandscape.size() - 1) {
			double x0 = bootLandscape[minIdx - 1].first, w0 = bootLandscape[minIdx - 1].second;
			double x1 = bootLandscape[minIdx].first,     w1 = bootLandscape[minIdx].second;
			double x2 = bootLandscape[minIdx + 1].first, w2 = bootLandscape[minIdx + 1].second;

			// Vertex of parabola through (x0,w0), (x1,w1), (x2,w2)
			double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
			double a = (x2 * (w1 - w0) + x1 * (w0 - w2) + x0 * (w2 - w1)) / denom;

			if (a > 0) {
				double b = (x2 * x2 * (w0 - w1) + x1 * x1 * (w2 - w0) + x0 * x0 * (w1 - w2)) / denom;
				double xStar = -b / (2.0 * a);

				// Clamp to the bracket to avoid extrapolation
				xStar = std::max(x0, std::min(x2, xStar));

				// Linearly interpolate IOD at xStar
				double iod0 = iodLandscape[minIdx - 1].second;
				double iod1_val = iodLandscape[minIdx].second;
				double iod2 = iodLandscape[minIdx + 1].second;
				if (xStar <= x1) {
					double t = (x1 - x0 > 0) ? (xStar - x0) / (x1 - x0) : 0.0;
					bootIOD = iod0 + t * (iod1_val - iod0);
				} else {
					double t = (x2 - x1 > 0) ? (xStar - x1) / (x2 - x1) : 0.0;
					bootIOD = iod1_val + t * (iod2 - iod1_val);
				}
			} else {
				bootIOD = iodLandscape[minIdx].second;
			}
		} else {
			bootIOD = iodLandscape[minIdx].second;
		}

		bootstrapIODs.push_back(bootIOD);
	}

	std::sort(bootstrapIODs.begin(), bootstrapIODs.end());
	double ciLow = bootstrapIODs[static_cast<int>(std::floor(0.025 * nBootstrap))];
	double ciHigh = bootstrapIODs[static_cast<int>(std::ceil(0.975 * nBootstrap)) - 1];

	std::cerr << "Parameter search complete." << std::endl << std::endl;
	std::cerr << "Summary:" << std::endl;
	std::cerr << "--------- " << std::endl;
	std::cerr << "Optimal fr = " << std::scientific << std::setprecision(3) << bestFr << std::endl;
	std::cerr << "Wasserstein distance = " << std::fixed << std::setprecision(3) << bestW << std::endl;
	std::cerr << "Median IOD = " << std::fixed << std::setprecision(0) << bestIOD << " kb" << std::endl;
	std::cerr << "95% confidence interval: [" << std::fixed << std::setprecision(0) << ciLow << " kb , " << std::fixed << std::setprecision(0) << ciHigh << " kb]" << std::endl;

	// Write output file
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto startTime = oss.str();
	std::ofstream outFile(args.output);
	outFile << "#DetectFile " << args.DetectInput << "\n";
	outFile << "#ForkFiles " << args.lForkInput << " " << args.rForkInput << "\n";
	outFile << "#OriginFile " << args.originInput << "\n";
	outFile << "#TerminationFile " << args.terminationInput << "\n";
	outFile << "#SystemStartTime " + startTime + "\n";
	outFile << "#Software " << std::string(getExePath()) << "\n";
	outFile << "#Version " << std::string(VERSION) << "\n";
	outFile << "#Commit " << std::string(getGitCommit()) << "\n";
	outFile << "#nForks " << data_behindDists.size() << "\n";
	outFile << "#MeanForkSpeed " << std::fixed << std::setprecision(3) << avgForkSpeed << " kb/min\n";
	outFile << "#OptimalFiringRate " << std::scientific << std::setprecision(6) << bestFr << "\n";
	outFile << "#WassersteinDistance " << std::fixed << std::setprecision(6) << bestW << "\n";
	outFile << "#MedianIOD " << std::fixed << std::setprecision(1) << bestIOD << " kb\n";
	outFile << "#95ConfidenceInterval " << std::fixed << std::setprecision(1) << ciLow << " kb , " << ciHigh << " kb\n";
	outFile << ">DataBehindDistances:\n";
	for (const auto& val : data_behindDists) {
		if (val >= 2.0) outFile << val << "\n";
	}
	outFile << ">SimBehindDistances:\n";
	for (const auto& val : bestSimDist.first) {
		if (val >= 2.0) outFile << val << "\n";
	}
	outFile << ">DataAheadDistances:\n";
	for (const auto& val : data_aheadDists) {
		if (val >= 2.0) outFile << val << "\n";
	}
	outFile << ">SimAheadDistances:\n";
	for (const auto& val : bestSimDist.second) {
		if (val >= 2.0) outFile << val << "\n";
	}

	// Write landscape from points collected during the search (sorted by fr)
	std::sort(landscapePoints.begin(), landscapePoints.end());
	outFile << ">Landscape:\n";
	for (const auto& pt : landscapePoints) {
		outFile << std::scientific << std::setprecision(6) << std::get<0>(pt) << "\t"
				<< std::fixed << std::setprecision(1) << std::get<1>(pt) << "\t"
				<< std::fixed << std::setprecision(4) << std::get<2>(pt) << "\n";
	}
	outFile.close();

	return 0;
}
