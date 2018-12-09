#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <unordered_set>
#include <thread>
#include <chrono>
#include <sys/stat.h>
using namespace std;

// compile with g++ -std=c++11 find_pauses.cpp 



// Function to split bedGraph into single basepair intervals 
int splitBedGraph( std::string bedGraphName, std::string bedGraphSplit ) {

    std::ifstream bedGraph( bedGraphName );
    ofstream outputFile( bedGraphSplit.c_str() );
    
    std::string line;
    std::string delim = "\t";
    size_t      pos   = 0;
 
    while (std::getline(bedGraph, line)) {
        std::vector<std::string> cols;
        std::string lineCol = line;

        // Create vector containing all columns for the line
        while ((pos = lineCol.find(delim)) != std::string::npos) {
            std::string col = lineCol.substr( 0, pos );
            cols.push_back( col );
            lineCol.erase( 0, pos + delim.length() );
        }

        cols.push_back( lineCol );  // add last column to the vector

        // Assign variables for each column
        std::string chrom = cols[0];
        int         start = std::atoi( cols[1].c_str() );
        int         end   = std::atoi( cols[2].c_str() );
        std::string count = cols[3];

        int length = end - start;

        // Divide bedGraph into single bp intervals 
        for (int i = 0; i < length; i++) {
            int startNew = start + i;
            int endNew   = startNew + 1;
        
            stringstream startStream;
            startStream << startNew;
            std::string startOut = startStream.str();              

            stringstream endStream;
            endStream << endNew;
            std::string endOut = endStream.str();

            std::string newLine = 
                chrom    + "\t" + 
                startOut + "\t" + 
                endOut   + "\t" + 
                count    + "\n";
            
            outputFile << newLine;
        }
    }

    bedGraph.close();

    return 0;
}





// Function to generate windows around each bedGraph interval 
int runBedtools( std::string bedGraphSplit, std::string geneRegion, std::string chromSizes, 
                 std::string leftStr, std::string rightStr, std::string bedGraphWins ) {

    std::string bashCommands = 
        "cat " + geneRegion + 
        " | sort -k1,1 -k2,2n" +
        " | bedtools intersect -sorted -a - -b " + bedGraphSplit + 
        " | awk -v OFS=\"\\t\" '{ $4 = $1\":\"$2\"-\"$3\"*\"$4; print }'" +
        " | bedtools slop -i - -g " + chromSizes + " -l " + leftStr + " -r " + rightStr + 
        " | sort -k1,1 -k2,2n -k3,3n -u" + 
        " | bedtools intersect -sorted -a - -b " + bedGraphSplit + " -wo" +
        " | sort -k4,4 -k1,1 -k2,2n -k3,3n" + 
        " > " + bedGraphWins;

    system( bashCommands.c_str() );

    return 0;  
}





// Function to identify pause sites based on countLim and stdevLim 
int findPauses( std::string bedGraphWins, int countLim, int stdevLim ) {

    // countLim is the min number of reads that must map to a bedGraph
    // interval for it to be considered

    // The signal must be (stdevLim * window stdev) greater than the
    // window mean for the interval to be called a pause site
 
    std::unordered_set<std::string> pauseSet;
    int pauseNum = -1; 

    while (pauseNum != 0) {
        pauseNum = 0;
        std::ifstream       Wins( bedGraphWins );
        std::string         oldName    = "";
        double              countTot   = 0;
        double              lenTot     = 0;
        int                 winSize    = 0;
        std::string         oldWinChrom;
        int                 pauseStart;
        int                 pauseEnd;
        int                 oldPauseStart;
        int                 oldPauseEnd;
        std::string         oldStrand;
        int                 pauseCount;
        std::vector<double> counts;
        std::string         line;

        // Read windows line by line
        while (std::getline( Wins, line )) {
            std::vector<std::string> cols;
            std::string delim = "\t";
            size_t      pos   = 0;

            // Create vector containing all columns for the line
            while ((pos = line.find(delim)) != std::string::npos) {
                std::string col = line.substr( 0, pos );
                cols.push_back( col );
                line.erase( 0, pos + delim.length() );         
            }

            cols.push_back( line );  // add last column to vector
        
            // Asign variables for each column
            std::string winChrom    = cols[0];
            std::string winStartStr = cols[1];
            std::string winEndStr   = cols[2];
            int         winStart    = std::atoi( cols[1].c_str() );
            int         winEnd      = std::atoi( cols[2].c_str() );
            std::string name        = cols[3];
            std::string strand      = cols[5];
            std::string bedStartStr = cols[7];
            std::string bedEndStr   = cols[8];
            int         bedStart    = std::atoi( cols[7].c_str() );
            int         bedEnd      = std::atoi( cols[8].c_str() );
            double      count       = std::atof( cols[9].c_str() ); 
            double      bedSize     = std::atof( cols[10].c_str() ); 

            winSize = winEnd - winStart;
         
            // Identify pause position within window
            pauseStart = winStart + ( winSize / 2 );
            pauseEnd   = pauseStart + 1;

            std::stringstream startStream;
            startStream << pauseStart;
            std::string pauseStartStr = startStream.str();
        
            std::stringstream endStream;
            endStream << pauseEnd;
            std::string pauseEndStr = endStream.str();

            // Check for pause coords in pauseSet
            std::string pauseCheck = winChrom + ":" + pauseStartStr + "-" + pauseEndStr + strand;

            if (pauseSet.find( pauseCheck ) != pauseSet.end()) {
                continue;
            } 

            // Check for bed coords in pauseSet
            std::string bedCheck = winChrom + ":" + bedStartStr + "-" + bedEndStr + strand;

            if (pauseSet.find( bedCheck ) != pauseSet.end()) {
                continue;
            }

            // Add up counts that are present in the same window
            if (oldName == "" || name == oldName) {
                counts.push_back( count ); 
                countTot += count;
                lenTot += bedSize;
            
                if (pauseStart == bedStart) {
                    pauseCount = count;
                }
            }

            // When the window changes calculate stats for the previous window
            else {
                if (pauseCount > countLim) {
                    double zeroLen   = winSize - lenTot;
                    double meanCount = countTot / winSize;
                    double sumVar = 0;

                    // Add zeros for bedGraph intervals that lacked signal
                    for (int i = 0; i < zeroLen; i++) {
                        counts.push_back( 0 );
                    }
            
                    // Calculate the mean and standard deviation for each window
                    for (int i = 0; i < counts.size(); i++) {
                        double var = counts[ i ] - meanCount; 
                        sumVar += pow( var, 2 );
                    }

                    double meanVar = sumVar / ( winSize - 1 );
                    double stdev   = sqrt( meanVar );
                    double limit   = meanCount + ( stdev * stdevLim );

                    // Idenitify bedGraph intervals where the number of counts is
                    // greater than the stdev * stdevLim
                    if (pauseCount > limit) {

                        pauseNum += 1;

                        // Pause stats 
                        std::stringstream countStream;
                        countStream << pauseCount;
                        std::string countStr = countStream.str();
                    
                        std::stringstream totStream;                
                        totStream << countTot;
                        std::string totStr = totStream.str();

                        std::string pauseStats = countStr + "," + totStr;

                        // Write pause coords to stdout
                        std::cout << oldWinChrom   << "\t" 
                                  << oldPauseStart << "\t" 
                                  << oldPauseEnd   << "\t" 
                                  << oldName       << "\t" 
                                  << pauseStats    << "\t"
                                  << oldStrand     << endl;

                        // Add pause coords to pauseSet
                        std::stringstream startStream;
                        startStream << oldPauseStart;
                        std::string pauseStartStr = startStream.str();
        
                        std::stringstream endStream;
                        endStream << oldPauseEnd;
                        std::string pauseEndStr = endStream.str();

                        std::string pauseCoords = oldWinChrom + ":" + pauseStartStr + "-" + pauseEndStr + oldStrand;
                        pauseSet.insert( pauseCoords );
                    }
                }

                // Reset everything for the next window
                counts.clear();
                countTot = 0;
                lenTot   = 0;

                // Begin operating on the new window
                counts.push_back( count ); 
                countTot += count;
                lenTot += bedSize;
            
                if (pauseStart == bedStart) {
                    pauseCount = count;
                }
            }

            oldName       = name;
            oldWinChrom   = winChrom;
            oldPauseStart = pauseStart;
            oldPauseEnd   = pauseEnd;
            oldStrand     = strand;
        }

        Wins.close();
    }

    return 0;
}





// Function to generate string used for naming output files 
std::string getName( std::string path ) {

    std::vector<string> pathVec;
    std::string pathDelim = "/";
    size_t pos = 0;
    
    while ((pos = path.find(pathDelim)) != std::string::npos) {
        std::string name = path.substr( 0, pos );
        pathVec.push_back( name );
        path.erase( 0, pos + pathDelim.length() );         
    }

    pathVec.push_back( path );  // add last filename to vector

    std::string lastName = pathVec[ pathVec.size() - 1 ];
    std::vector<string> nameVec;
    std::string nameDelim = ".";

    while ((pos = lastName.find(nameDelim)) != std::string::npos) {
        std::string name = lastName.substr( 0, pos );
        nameVec.push_back( name );
        lastName.erase( 0, pos + nameDelim.length() );
    }

    std::string sampleName;

    for (int i = 0; i < nameVec.size(); i++) {
        sampleName += nameVec[ i ];        

        if (i != nameVec.size() - 1) {
            sampleName += ".";
        }
    }
    
    return sampleName;
}





// Main function
int main( int argc, char *argv[] ) {

    // Input arguments
    std::string bedGraph   = argv[1];
    std::string geneRegion = argv[2];
    std::string chromSizes = argv[3];
    std::string countLimIn = argv[4];
    std::string stdevLimIn = argv[5];
    int countLim           = std::atoi( countLimIn.c_str() );
    int stdevLim           = std::atoi( stdevLimIn.c_str() );

    // Extract bedGraph name  
    std::string bedGraphName  = getName( bedGraph );
    
    // Temporary files  
    std::string bedGraphSplit = bedGraphName + ".split";    
    std::string bedGraphWins  = bedGraphName + ".wins";

    // Split bedGraph into single base pair intervals 
    splitBedGraph( bedGraph, bedGraphSplit );

    // Generate windows around each bedGraph interval  
    runBedtools( bedGraphSplit, geneRegion, chromSizes, "100", "99", bedGraphWins );

    // Identify pause sites 
    findPauses( bedGraphWins, countLim, stdevLim );

    // Remove temporary files
    std::string rmCommands = "rm " + bedGraphSplit + " " + bedGraphWins; 
    system( rmCommands.c_str() ); 

    return 0;
}





