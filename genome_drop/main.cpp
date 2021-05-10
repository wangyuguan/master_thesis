/* file: main.cpp
 * project: GenomeDrop
 * copyright Mark Abney, 2013
 *
 * This program drops a genome, consisting of a specified number of chromosomes,
 * down a known pedigree. The genotypes of a study sample at regularly spaced
 * markers separated by a given recombination fraction are output to a file.
 * Optionally, IBD matrices for each marker and their Cholesky decompositions are
 * also output into a binary formatted file.
 */

#define VERSION 0.1

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>

#include "genedrop.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    int opt;
    string basename;
    string paramfile;
    string outfile;
    string logfile;
    string ibdfile;
    string seedstr;
    bool output_ibd = false;
    long seed = 0;

    while ((opt = getopt(argc, argv, ":p:io:s:")) != -1) {
        switch (opt) {
            case 'p':
                paramfile = optarg;
                break;
            case 'i':
                output_ibd = true;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 's':
                seedstr = optarg;
                seed = strtol(optarg, (char **) NULL, 10);
                break;
            case ':':
                cout << "Option " << (char)optopt << " requires an argument." << endl;
                exit(EXIT_FAILURE);
            case '?':
                cout << "Unknown option: " << (char)optopt << endl;
                exit(EXIT_FAILURE);
            default:
                cerr << "Unexpected error.\n";
                exit(EXIT_FAILURE);
        }
    }
       
    if (paramfile.empty()) {
        cerr << "A parameter file must be specified.\n";
        exit(EXIT_FAILURE);
    } else {
        cout << "parameter file: " << paramfile << endl;
    }

    if (seedstr.empty()) {
        cerr << "An integer seed value must be specified.\n";
        exit(EXIT_FAILURE);
    }

    if (outfile.empty()) {
        string::size_type idx = paramfile.rfind('.');
        if (idx == string::npos || idx == 0) {
            basename = paramfile;
        } else {
            basename = paramfile.substr(0, idx);
        }
        basename += "_" + seedstr;
    } else {
        basename = outfile;
    }
    outfile = basename + ".geno";
    //logfile = basename + '_' + seedstr + ".log";
    logfile = basename + ".log";

    if (output_ibd)
        ibdfile = basename + "." + "ibdm";
    cout << "output file: " << outfile << endl;
    cout << "log file: " << logfile << endl;

    ofstream logfl(logfile.c_str());
    if ( !logfl) {
        cerr << "Cannot open log file: " << logfile << endl;
        exit(EXIT_FAILURE);
    }
    logfl << "GenomeDrop Version: " << VERSION << endl;
    for (int i = 0; i < argc; ++i) {
        logfl << argv[i] << " ";
    }
    logfl << endl;
    logfl << "Random number seed: " << seed << endl;
    logfl << "Input parameter file: " << paramfile << endl;
    logfl << "Output genotype file: " << outfile << endl;
    if (output_ibd)
        logfl << "Output IBD matrices file: " << ibdfile << endl;
    else
        logfl << "No output IBD matrices file." << endl;

    /*
     * Parse the parameter file and set up necessary data structures, etc.
     */
    GeneDrop genomedrop(paramfile, outfile, output_ibd, logfl);
    genomedrop.setup();
    genomedrop.random_init(seed);
    genomedrop.simulate();
    genomedrop.print_genome();
    if (output_ibd)
        genomedrop.print_ibd(ibdfile);

}
