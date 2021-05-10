/* file: genedrop.cpp
 * Project: GenomeDrop
 * copyright Mark Abney, 2013
 *
 * The necessary classes and codes for doing a gene dropping simulation.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>
#include <random>

#include "genedrop.hpp"

using namespace std;

void log_mesg(ostream &logstm, const string &msg) {
    logstm << msg << endl;
    cerr << msg << endl;
}

template <typename T>
T get_next_param(ifstream &inputfl, const string &msg, ostream &outfl=cerr)//, string &output) {
{
    T output;
    string line;
    string::size_type pos;
    do {
        getline(inputfl, line);
    	pos = line.find_first_not_of(" \t");
    } while (line.empty() || pos == string::npos || line[pos] == '#') ;
    //getline(inputfl, line);
    stringstream(line) >> output;
    //outfl << msg << output << endl;
    ostringstream omsg;
    omsg << msg << output;
    log_mesg(outfl, omsg.str());
    //cout << msg << output << endl;
    return output;
}


// int Pedigree::get_id( string &findiv) const {
//         IdMap::const_iterator it = id_of.find(findiv);
//         if (it != id_of.end()) {
//             return it->second;
//         } else {
//             return -1;
//         }
// }

void Pedigree::readped(const string &pedfile, ofstream &logfl) {
    /* The pedigree file should have columns: individual, mother, father, sex, ...
     * Parents must come before children and non-founders must have both parents.
     * Founders should have the character '0' as both parents.
     * Lines that are empty or begin with # are skipped.
     */
    cout << "Reading " << pedfile << endl;
    ifstream pedfl(pedfile.c_str());
    if (! pedfl) {
        log_mesg(logfl, "Cannot open pedigree file: " + pedfile);
        //cerr << "Cannot open pedigree file: " << pedfile << endl;
        //logfl << "Cannot open pedigree file: " << pedfile << endl;
        exit(10);
    }

    bool pederr = false;
    id_of["0"] = 0;
    findiv_of.push_back("0");
    isfounder.push_back(true);
    // DEBUG
    // cout << "findiv_of[0] = " << findiv_of[0] << endl;
    // cout << "id_of 0 is: " << id_of["0"] << ";\n";
    // cout << "size of parent = " << parent.size() << ", capacity = " << parent.capacity() << endl;
    // END DEBUG
    array<int,2> init_ar = {0,0};
    parent.push_back(init_ar);
    // DEBUG
    // cout << "at 0,0: " << parent.at(0)[0] << endl;
    // cout << "size of parent = " << parent.size() << ", capacity = " << parent.capacity() << endl;
    // END DEBUG
    //parent[0][0]  = 0; //Founder parents also have the founder ID.
    string line;
    for (int id = 0; getline(pedfl, line); ) {
        //array<int,2> x = {0,0};
        //parent.push_back(x);
	string::size_type pos = line.find_first_not_of(" \t");
        if (line.empty() || pos == string::npos || line[pos] == '#') {
            //if (line.empty() || line[0] == '#') {
            cout << "skipping: " << line << endl;
            continue;
        }
        // Debug:
        // cout << line << endl;
        // End debug
        ++npeople;
        ++id;
        parent.push_back(init_ar);
        // DEBUG
        // cout << "size of parent = " << parent.size() << ", capacity = " << parent.capacity() << endl;
        // END DEBUG
        
        //int mid=0, fid=0;
        string findiv, mother, father, sex;
        stringstream(line) >> findiv >> mother >> father >> sex;
        // DEBUG
        // cout <<"id :" << id << ": findiv :" << findiv << ": mother :" << mother << ": father :" << father
        // << ": sex :" << sex <<":"<<endl;
        // END DEBUG
        if (id_of.count(findiv) > 0) {
            log_mesg(logfl, "Error: Person " + findiv + " appears multiple times in the pedigree.");
            //cerr << "Error: Person " << findiv << " appears multiple times in the pedigree." << endl;
            //logfl << "Error: Person " << findiv << " appears multiple times in the pedigree." << endl;
            pederr = true;
        }
        if (id_of.count(mother) == 0) {
            log_mesg(logfl, "Error: Person " + findiv + " appears before parent " + mother);
            //cerr << "Error: Person " << findiv << " appears before parent " << mother << endl;
            //logfl << "Error: Person " << findiv << " appears before parent " << mother << endl;
            pederr = true;
        } else {
            //mid = id_of[mother];
            //cout << "Adding mother: :" << mother <<":" << endl;
            parent[id][0] = id_of[mother];
	}
        if (id_of.count(father) == 0) {
            log_mesg(logfl, "Error: Person " + findiv + " appears before parent " + father);
            //cerr << "Error: Person " << findiv << " appears before parent " << father << endl;
            //logfl << "Error: Person " << findiv << " appears before parent " << father << endl;
            pederr = true;
        } else {
            //fid = id_of[father];
            //cout << "Adding father: :" << father << ":" << endl;
            parent[id][1] = id_of[father];
	}
        // Individuals are not allowed to have only one unknown parent
        if ((mother=="0" && father!="0") || (mother!="0" && father=="0")) {
            log_mesg(logfl, "Error: Person " + findiv + " has only one unknown parent.");
            //cerr << "Error: Person " << findiv << " has only one unknown parent." << endl;
            //logfl << "Error: Person " << findiv << " has only one unknown parent." << endl;
            pederr = true;
        } /*else if ( ! mother && ! father ) {
            ++nfounders;
        } else {
            nonfounder[nf++] = id;
            }*/
        
        // cout << "Adding findiv: :" << findiv << ": id is :" << id << ":" << endl;
        id_of[findiv] = id;
	//findiv_of[id] = findiv;
        findiv_of.push_back(findiv);
        if (mother=="0" && father=="0")
        {
            isfounder.push_back(true);
            ++nfounder;
        } else {
            isfounder.push_back(false);
        }
	// Debug:
        // cout << "id_of [findiv] and findiv_of[id]" << endl;
	// cout << id_of[findiv] << "(" << findiv_of[id] << ") "<< parent[id][0];
	// cout << "(" << findiv_of[parent[id][0]] << ") ";
	// cout << parent[id][1] << "(" << findiv_of[parent[id][1]] <<")" << endl;
	// End Debug
    }
    pedfl.close();
    
    if (pederr) {
        log_mesg(logfl, "Exiting due to pedigree errors.");
        //cerr << "Exiting due to pedigree errors." << endl;
        //logfl << "Exiting due to pedigree errors." << endl;
        exit(11);
    }
    log_mesg(logfl, "Successfully read pedigree.");
}


int HaploData::n_ibd(const int chr, const int loc, const int id1, const int id2) {
    int count = 0;
    for (int p1 = 0; p1 < 2; ++p1) {
        for (int p2 = 0; p2 < 2; ++p2) {
            if (haplo[id1][chr][loc][p1] == haplo[id2][chr][loc][p2])
                count++;
        }
    }
    return count;
}


void GeneDrop::print_ibd(const string &ibdfile) {
    ostream* outstrp = &cout;
    ofstream ibdfl(ibdfile.c_str());
    if (! ibdfl) {
        log_mesg(logfl, "Cannot open ibd output file: " + ibdfile + "\nUsing stdout....");
    } else {
        outstrp = &ibdfl;
    }
    ostream& ibdout = *outstrp;

    ibdout << "#chrom\tsnp\tcM\tbp";
    for (int idx1 = 0; idx1 < study_sample.size(); ++idx1) {
        for (int idx2 = 0; idx2 <= idx1; ++idx2) {
            ibdout << "\t" << pedigree.get_findiv(study_sample[idx1]) << "," << pedigree.get_findiv(study_sample[idx2]);
        }
    }
    ibdout << endl;

    double bp_pos = 0; // Put this here to give base pair position from beginning of genome
    for (int chr = 0; chr < genome.n_chromosomes(); ++chr) {
        double cMspace = genome.cM_space(chr);
        string chrname = genome.get_name(chr);
        for (int loc = 0; loc < genome.n_markers(chr); ++loc) {
            ibdout << chrname << "\tc" << chrname << "s" << loc << "\t" << loc * cMspace << "\t" << bp_pos;
            bp_pos += cMspace * 1e6; // Assumes 1 cM = 1 Mb
            for (int idx1 = 0; idx1 < study_sample.size(); ++idx1) {
                int id1 = study_sample[idx1];
                for (int idx2 = 0; idx2 <= idx1; ++idx2) {
                    int id2 = study_sample[idx2];
                    int n_ibd = haplo.n_ibd(chr, loc, id1, id2);
                    ibdout << "\t" << 0.25 * n_ibd;
                }
            }
            ibdout << endl;
        }
    }

    if (ibdfl)
        ibdfl.close();
}

void GeneDrop::print_genome() {
    ostream* true_outp = &cout;
    ostream* geno_outp = &cerr;

    ofstream genofl(outfile.c_str());
    if (! genofl) {
        log_mesg(logfl, "Cannot open genotype file: " + outfile + "\nUsing stderr....");
    } else {
        geno_outp = &genofl;
    }
    ostream& genoout = *geno_outp;
    
    string truename = outfile + "_true"; // This file will have genotypes of founder alleles
    ofstream truefl(truename.c_str());
    if ( ! truefl ) {
        log_mesg(logfl, "Cannot open output file: " + truename + "\nUsing stdout....");
        //truefl = cout;
    } else {
        true_outp = &truefl;
    }
    ostream& trueout = *true_outp;

    trueout << "#chrom\tsnp\tcM\tbp";
    genoout << "#chrom\tsnp\tcM\tbp";
    for (int sid = 0; sid < study_sample.size(); ++sid) {
        genoout << "\t" << pedigree.get_findiv(study_sample[sid]);
        trueout << "\t" << pedigree.get_findiv(study_sample[sid]);
    }
    genoout << endl;
    trueout << endl;
    double bp_pos = 0; // Put this here to give base pair position from beginning of genome
    for (int chr = 0; chr < genome.n_chromosomes(); ++chr) {
        double cMspace = genome.cM_space(chr);
        string chrname = genome.get_name(chr);
        for (int loc = 0; loc < genome.n_markers(chr); ++loc) {
            haplo.assign_types(chr, loc, rand_engine);
            genoout << chrname << "\tc" << chrname << "s" << loc << "\t" << loc * cMspace << "\t" << bp_pos;
            trueout << chrname << "\tc" << chrname << "s" << loc << "\t" << loc * cMspace << "\t" << bp_pos;
            bp_pos += cMspace * 1e6; // Assumes 1 cM = 1 Mb
            for (int sid = 0; sid < study_sample.size(); ++sid) {
                int id = study_sample[sid];
                int ale1, ale2;
                trueout << "\t" << (ale1 = haplo.get_allele(id, chr, loc, 0))
                        << " " << (ale2 = haplo.get_allele(id, chr, loc, 1));
                genoout << "\t" << haplo.allelic_type_of(chr, loc, ale1) << " "
                        << haplo.allelic_type_of(chr, loc, ale2);
            }
            genoout << endl;
            trueout << endl;
        }
    }
    if (genofl)
        genofl.close();
    log_mesg(logfl, "Successfully output the genotypes.");
    if (truefl)
        truefl.close();
    log_mesg(logfl, "Successfully output the true genotypes.");
}

void GeneDrop::random_init(int rseed = 0) {
    unsigned defseed = 5489; // Proper default seed for mt19937
    rand_engine.seed(defseed + rseed);
}

void GeneDrop::setup() {
    // Just parse the parameter file and get needed variables.

    
    ifstream inputfl(paramfile.c_str());
    if ( ! inputfl ) {
        log_mesg(logfl, "Cannot open parameter file: " + paramfile + "\nExiting...");
        //cerr << "Cannot open parameter file: " << paramfile << endl;
        //logfl << "Cannot open parameter file: " << paramfile << endl;
        //logfl << "Exiting..." << endl;
        exit(EXIT_FAILURE);
    }

    string pedfile = get_next_param<string>(inputfl, "Pedigree file: ", logfl);
    string studyfile = get_next_param<string>(inputfl, "Study sample file: ", logfl);
    string mapfile = get_next_param<string>(inputfl, "Map file: ", logfl);
    //mafreq = get_next_param<double>(inputfl, "Minor allele freq: ", logfl);
    miss_rate = get_next_param<double>(inputfl, "Missing genotype rate: ", logfl);
    allele_err = get_next_param<double>(inputfl, "Allelic error rate: ", logfl);
    inputfl.close();
    
    pedigree.readped(pedfile, logfl);
    // Debug:
    // for (int i = 0; i < pedigree.n_people(); ++i) {
    // cout << "id: " << i << " findiv: " << pedigree.get_findiv(i) << " founder: " << pedigree.is_founder(i) << endl;
    // }
    // End debug
    readsample(studyfile);
    readmap(mapfile);
    // DEBUG
    // cout << "First address of genome: " << &genome << endl;
    // cout << "First address of pedigree: " << &pedigree << endl;
    // END DBUG
    haplo.init(genome, pedigree);
}

void GeneDrop::simulate() {
    // Here we do the gene dropping simulation. The plan is to loop through every person in
    // the pedigree and assign a value for each marker on each of the person's haplotypes. If
    // the person is a founder, the two haplotypes have constant values for every marker this
    // identifies each haplotype uniquely. Fon non-founders the maternal haplotype gets one of the
    // values from the person's mother based on Mendelian inheritance and recombination probabilities.
    int npeople = pedigree.n_people();
    int founder_count = 0;
    bernoulli_distribution bern_val(0.5);
    uniform_real_distribution<double> rand_unif(0, 1);
    for (int person = 1; person <= npeople; ++person) {
        if (pedigree.is_founder(person)) {
            ++founder_count;
            for (int chr=0; chr < genome.n_chromosomes(); ++chr) {
                for (int loc = 0; loc < genome.n_markers(chr); ++loc) {
                    haplo.set_allele(person, chr, loc, 0) = 2 * founder_count - 1;
                    haplo.set_allele(person, chr, loc, 1) = 2 * founder_count;
                }
            }
        } else { // person is not a founder:
            for (int chr=0; chr < genome.n_chromosomes(); ++chr) {
                int nloci = genome.n_markers(chr);
                for (int parental_allele = 0; parental_allele < 2; ++parental_allele) {
                    // For each parental allele of the person decide which grandparental chromosome is transmitted
                    // i.e. which of the two parental chromosomes is transmitted (= trchrom)
                    bool trchrom = bern_val(rand_engine);
                    int parent = pedigree.parent_of(person, parental_allele);
                    haplo.set_allele(person, chr, 0, parental_allele) = haplo.get_allele(parent, chr, 0, trchrom);
                    for (int locus = 1 ; locus < nloci; ++locus) {
                        // And now proceed through all the markers allowing for possible recombination
                        trchrom = rand_unif(rand_engine) < genome.recomb_frac(chr) ? 1 - trchrom : trchrom;
                        haplo.set_allele(person, chr, locus, parental_allele) = haplo.get_allele(parent, chr, locus, trchrom);
                    } // end marker loop
                } // end parent loop
            } // end chromosome loop
        } // end person founder or not condition
    } // end person loop
    log_mesg(logfl, "Completed gene dropping simulation.");
}

void GeneDrop::readmap(const string &mapfile) {
    // Read the map file. Rows are chromosome and columns are chromosome name,
    // length (in cM), marker spacing (in cM), allele frequency for the first allele
    // of a biallelic marker.
    ifstream infl(mapfile.c_str());
    if (! infl) {
        log_mesg(logfl, "ERROR: Cannot open map file: " + mapfile);
        exit(13);
    }

    log_mesg(logfl, "Reading map file: " + mapfile);
    bool maperr = false;
    string line;
    for (int chr_id = 0; getline(infl, line); ) {
	string::size_type pos = line.find_first_not_of(" \t");
        if (line.empty() || pos == string::npos || line[pos] == '#') {
            cout << "skipping: " << line << endl;
            continue;
        }
	// Read each chromosome and build up the Genome object:
        string name;
        double len, spac, maf;
        stringstream(line) >> name >> len >> spac >> maf;
        if (len < 0) {
            log_mesg(logfl, "ERROR: Chromosome length " + to_string(len) + " for chromosome "
                     + name + " is invalid.");
            exit(14);
        }
        if (spac <= 0) {
            log_mesg(logfl, "ERROR: marker spacing " + to_string(spac) + " for chromosome "
                     + name + " is invalid.");
            exit(15);
        }
        if (maf <= 0) {
            log_mesg(logfl, "ERROR: Allele frequency " + to_string(maf) + " for chromosome "
                     + name + " is invalid.");
            exit(16);
        }

        genome.add_chromosome(name, len, spac, maf);
    }
    infl.close();
    log_mesg(logfl, "Genome has " + to_string(genome.n_chromosomes()) + " chromosomes.");
}

void GeneDrop::readsample(const string &studyfile) {
    // Read the study sample as a list of findivs that are white space delimited.
    ifstream infl(studyfile.c_str());
    if (! infl) {
        log_mesg(logfl, "ERROR: Cannot open study sample file: " + studyfile);
        exit(10);
    }

    bool studyerr = false;
    string findiv;
    while (infl >> findiv) {
        int id = pedigree.get_id(findiv);
        if (id < 0) {
            log_mesg(logfl, "ERROR: Study sample ID " + findiv + " is not in the pedigree.");
            studyerr = true;
        } else {
            study_sample.push_back(id);
        }
    }
    infl.close();
    
    if (studyerr) {
        log_mesg(logfl, "Errors in the study sample.\nExiting....");
        exit(12);
    }
    log_mesg(logfl, "Number in the study sample: " + to_string(study_sample.size()));
    //cerr << "Number in the study sample: " << study_sample.size() << endl;
    //logfl << "Number in the study sample: " << study_sample.size() << endl;
    // Debug:
    // for (int i=0; i < study_sample.size(); ++i)
    // cout << i << " " << study_sample[i] << " " << pedigree.get_findiv(study_sample[i]) << endl;
    // End Debug
}
