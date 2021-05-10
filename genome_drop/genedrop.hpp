/* file: genedrop .h
 * Project GenomeDrop
 * copyright Mark Abney, 2013
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>
#include <random>

//using namespace std;

class Chromosome {
  private:
    std::string name;
    double length; // in centiMorgans
    double spacing; // distance between markers in centiMorgans
    unsigned nmarkers;
    double allele_freq;
    double recomb_fraction;
    
  public:
    Chromosome() {}
    Chromosome(const std::string &nam, const double len, const double spac, const double maf)
            : name(nam), length(len), spacing(spac), allele_freq(maf)
    {
        //nmarkers = static_cast<unsigned> ( floor(length / spacing) + 1 );
        nmarkers = (unsigned) (floor(length / spacing) + 1);
        //cout << "Number of markers: " << n_markers << endl;
        recomb_fraction = 0.5 * tanh(2 * spacing); // Kosambi map function in case of a large spacing
    }
    unsigned n_markers() const { return nmarkers; }
    double recomb_frac() const { return recomb_fraction;}
    double cM_space() const {return spacing;}
    const std::string& get_name() const { return name; }
    double minor_allele_freq() const { return allele_freq; }
};

class Genome {
  public:
    void add_chromosome(const std::string &name, const double len, const double spac, const double maf) {
        chrom.push_back(Chromosome(name, len, spac, maf));
    }
    int n_chromosomes() const { return chrom.size(); }
    int n_markers(const int chr) const { return chrom[chr].n_markers(); }
    double recomb_frac(const int chr) const { return chrom[chr].recomb_frac(); }
    double cM_space(const int chr) const { return chrom[chr].cM_space(); }
    const std::string& get_name(const int chr) const { return chrom[chr].get_name(); }
    double minor_allele_freq(const int chr) const { return chrom[chr].minor_allele_freq(); }
    
  private:
    std::vector<Chromosome> chrom;
};

class Pedigree {
    typedef std::unordered_map<std::string, int> IdMap;
    
  public:
    Pedigree() : npeople(0), nfounder(0) {}
    void readped(const std::string &pedfile, std::ofstream &logfl);
    int get_id(const std::string &findiv) const {
        IdMap::const_iterator it = id_of.find(findiv);
        if (it != id_of.end()) {
            return (*it).second;
        } else {
            return -1;
        }
    }
    
    std::string get_findiv(const int id) const {
        if (id >= findiv_of.size() || id < 0)
            return "";
        else
            return findiv_of[id];
    }
    
    int n_people() const { return npeople;}
    bool is_founder(const int& id) const {
        return isfounder[id];
    }
    int parent_of(const int id, const int par) const {
        // par is either 0 for the mother or 1 for the father
        return parent[id][par];
    }
    int n_founders() const { return nfounder; }
    
  private:
    IdMap id_of; // The value is the internal ID for the key (findiv, i.e. file ID).
    std::vector<std::string> findiv_of; // Index is the internal ID, value is the findiv.
    std::vector<std::array<int,2>> parent; //parent[i][0] is mother of internal ID i, parent[i][1] is the father
    std::vector<bool> isfounder; //Indicator of whether ID is a founder
    int npeople;
    int nfounder;
};

class HaploData {
  public:
    HaploData() {}
    //HaploData(Genome &gen, Pedigree &ped) : genome(gen), pedigree(ped) {
    void init(const Genome &gen, const Pedigree &ped) {
        pgenome = &gen;
        ppedigree = &ped;
        const Genome &genome = *pgenome;
        const Pedigree &pedigree = *ppedigree;

        // DEBUG
        // std::cout << "Address of genome: " << &genome << std::endl;
        // std::cout << "Address of pedigree: " << &pedigree << std::endl;
        // END DEBUG
        
        int npeop = pedigree.n_people();
        int nchrom = genome.n_chromosomes();
        haplo.resize(npeop + 1); // Add 1 because IDs range from 1 to npeop (ID 0 is special)
        alleletype.resize(nchrom);
        for (int id = 0; id < npeop + 1; ++id) {
            haplo[id].resize(nchrom);
            for (int chr = 0; chr < nchrom; ++chr) {
                int nmark = genome.n_markers(chr);
                haplo[id][chr].resize(nmark);
                alleletype[chr].resize(nmark);
                for (int marker = 0; marker < nmark; ++marker) {
                    haplo[id][chr][marker].resize(2);
                }
            }
        }
    }

    // set_allele and get_allele either set or return founder alleles.
    int &set_allele(const int id, const int chr, const int loc, const int parent) {
        return haplo[id][chr][loc][parent];
    }
    const int get_allele(const int id, const int chr, const int loc, const int parent) const {
        return haplo[id][chr][loc][parent];
    }
    // allelic_type_of returns the observed allele (assuming no error) for the given founder allele
    const int allelic_type_of(const int chr, const int loc, const int f_ale) const {
        return alleletype[chr][loc][f_ale];
    }
    template <typename RandEng>
    void assign_types(const int chr, const int loc, RandEng& reng);

    int n_ibd(const int chr, const int loc, const int id1, const int id2);
    
  private:
    const Genome *pgenome;
    const Pedigree *ppedigree;
    std::vector<std::vector<std::vector<std::vector<int> > > > haplo; // haplo[ID][chrom][locus][maternal(0) or paternal(1)]
    std::vector<std::vector<std::vector<int> > > alleletype; // Allelic type of each founder allele for each locus and chrom
    // alleletype[chrom][loc][founder_allele] evaluates to a 1 or 2 (2 = minor allele)
};

template <typename RandEng>
void HaploData::assign_types(const int chr, const int loc, RandEng& reng) {
    const Genome &genome = *pgenome;
    const Pedigree &pedigree = *ppedigree;
    std::uniform_real_distribution<double> rand_unif(0, 1);
    int nfounders = pedigree.n_founders();

    alleletype[chr][loc].resize(2 * nfounders + 1);
    alleletype[chr][loc][0] = 0;
    for (int fale = 1; fale <= 2 * nfounders; ++fale) {
        double frac = rand_unif(reng);
        if (frac <= genome.minor_allele_freq(chr)) {
            alleletype[chr][loc][fale] = 2;
        } else {
            alleletype[chr][loc][fale] = 1;
        }
    }
}

class GeneDrop {
  public:
    //GeneDrop() {}
    GeneDrop(std::string &pfname, std::string &ofname, bool oibd, std::ofstream &ostr)
            : paramfile(pfname), outfile(ofname), output_ibd(oibd), logfl(ostr)
    {}
    void setup();
    void random_init(int rseed);
    void readsample(const std::string &studyfile);
    void readmap(const std::string &mapfile);
    void simulate();
    void print_genome();
    void print_ibd(const std::string &ibdfile);
    
  private:
    std::string paramfile;
    std::string outfile;
    bool output_ibd;
    std::ofstream &logfl;
    Pedigree pedigree;
    std::vector<int> study_sample;
    Genome genome;
    HaploData haplo;
    std::mt19937 rand_engine;
    //double mafreq;
    double miss_rate;
    double allele_err;
};
