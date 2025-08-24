#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <random>
#include <iterator>
//#include <boost/math/distributions/hypergeometric.hpp>
//#include <boost/math/policies/policy.hpp>

using namespace std;

/*
// hypergeometric test
// r: number of success
// n: sample size
// N: Total object
// k: number of success events drawn from n samples

double hypergeometrictest(int r, int n, int N, int k)
{
    if (k <= 0 )
        return 1;
    double cdf;

    
    boost::math::hypergeometric_distribution<double> hg_dist(r, n, N);

    
    cdf = boost::math::cdf<double>(hg_dist, k-1);

    return 1-cdf;

} 
*/

double overlap_coefficient( int sa, int sb, int si )
{
    if ( sa == 0 || sb == 0 )
    {
        cout<<"error cal oc: "<<sa<<" "<<sb<<endl; exit(1); 
    }

    int mi = sa;
    if ( sb < sa )
        mi = sb;
    double oc = (1.0*si)/mi;
    return oc;
}

// shaffling
void shaffle_vector( int vesize, int samplesize, int shaffle_times, vector<vector<int > > &shaffled_index ) 
// vesize: total vector size;
// samplesize: sampling vector size
// shaffling times
//void shaffle_vector()
{
    random_device rd;
    
    
    srand(time(NULL));
    vector<int > starts;
    int maxstart = vesize - samplesize;
    for ( int k = 0; k < shaffle_times; ++k )
    {
        int s = rand() % maxstart;
        starts.push_back(s);
    }
   
    for ( int k = 0; k < shaffle_times; ++k )
    {   
       
        vector<int> myvector;
        for (int i=0; i<vesize; ++i) 
            myvector.push_back(i); 
        
        mt19937 g(rd());
        shuffle( myvector.begin(), myvector.end(), g );
        
        int s = starts[k];
        vector<int> tmp;
     //   cout<<s<<endl;
        for ( int i = s; i < s+samplesize; ++i )
        {
            tmp.push_back(myvector[i]);
        }
        shaffled_index.push_back(tmp);
    }
}

void func1()
{
    vector<vector<int > > shaffled_index;

    shaffle_vector( 100, 10, 19, shaffled_index);
    for ( size_t i = 0; i < shaffled_index.size(); ++i )
    {
        for ( size_t j = 0; j < shaffled_index[i].size(); ++j )
            cout<<shaffled_index[i][j]<<" ";
        cout<<endl;
    }
}

void get_sampled_geneset( vector<string > &gs, vector<int > &index, set<string > &sgs )
{
    
    for ( size_t i = 0; i < index.size(); ++i )
    {
        if ( (int)index[i] >= gs.size() ) 
        {
            cout<<"error index i == "<<index[i] <<" > "<<gs.size()<<endl; exit(1);
        }
        sgs.insert( gs[index[i]]);
    }
}

void cal_mc_p_larger( int obs, vector<int > &rd_sampled, double &p, double &fc )
{
    int l = 0;
    int s = 0;
    int agg = 0;
    for ( int i = 0; i < (int)rd_sampled.size(); ++i )
    {
        if ( rd_sampled[i] >= obs )
            l +=1;
        else
            s += 1;
        agg += rd_sampled[i];
    }
    double ave = agg*1.0/(int)rd_sampled.size();
    p = l*1.0/(l+s);
    fc = (obs+1.0) / (ave+1.0);

}

vector<string > parse_string( string & instr, char spl )
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == spl )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
            if ( s != "" || instr[i] != ' ' )
    			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

vector<string > parse_string( string & instr)
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == '\t' || instr[i] == ' ' )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

string inttostr(int i )
{
	string Res;
	ostringstream convert;
	convert << i;
	Res = convert.str();
	return Res;
}

void readingeneset( string infile, set<string > &gs )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );

        string g = ps[0];
        gs.insert(g);
    }
    inf.close();
}

void readin_lrp( string infile, map<string, set<string > > &lrs, map<string, set<string > > &rls )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string l = ps[0];
        string r = ps[1];
        lrs[l].insert(r);
        rls[r].insert(l);
    }
    inf.close();

}

void readin_lrp( string infile, 
    map<pair<string, string>, set<pair<string, string> > > &lrs, 
    map<pair<string, string>, set<pair<string, string> > > &rls )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string cta = ps[0];
        string l = ps[1];
        string ctb = ps[2];
        string r = ps[3];
        lrs[make_pair(cta,l)].insert(make_pair(ctb,r));
        rls[make_pair(ctb,r)].insert(make_pair(cta,l));
    }
    inf.close();

}

void readin_lrp2( string infile, map<string, set<pair<string, string > > > &ligand_sender_receiver )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf,line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string cta = ps[0];
        string l = ps[1];
        string ctb = ps[2];
        string r = ps[3];
        ligand_sender_receiver[l].insert( make_pair(cta, ctb ) );

    }
    inf.close();
}

class Links
{
public:
    string nodea;
    string nodeb;
    string celltypea;
    string celltypeb;
    double cor;
    double pvalue;
    set<size_t > CC_ind;   // may be redundant
    set<string > CC_cps;   // unique interactions
    set<string > overlapped_genes;

    bool sig_ovl_gene_inct;
    bool sig_ovl_gene_crct;
    bool shared_gokegg;
    bool CC;

    Links(string innodea, string innodeb, double incor, double inp )
    {
        nodea = innodea;
        nodeb = innodeb;
        cor = incor;
        pvalue = inp;
        sig_ovl_gene_inct = false;
        sig_ovl_gene_crct = false;
        shared_gokegg = false;
        CC = false;
    }
};

class CellChat
{
public:
    string id;
    string partnera;
    string partnerb;
    string genea;
    string geneb;
    string classification;
    string celltypea;
    string celltypeb;
    string site;
    string broad_celltypea;
    string broad_celltypeb;
    set<string > geneseta;
    set<string > genesetb;
    CellChat( string inid, string inpa, string inpb, string inga, string ingb,
        string incl, string inca, string incb, string insite, string inbca, string inbcb )
    {
        id = inid;
        partnera = inpa;
        partnerb = inpb;
        genea = inga;
        geneb = ingb;
        classification = incl;
        celltypea = inca;
        celltypeb = incb;
        site = insite;
        broad_celltypea = inbca;
        broad_celltypeb = inbcb;
    }
    CellChat( string ingenea, string ingeneb )
    {
        genea = ingenea;
        geneb = ingeneb;
    }
};

class Bank
{
public:
    map<string, set<string > > cons_geneset;
    map<string, set<string > > module_geneset;

    void readinconsgeneset(string infile );
    void readinmodulegeneset( string infile );

    set<string> getgeneset( string item, string fea  );
    
};

class GeneSet
{
public:
    map<string, string > items_celltype;   // this is broad cell type
    map<string, string > items_feature;

    vector<Links> link_ve;
    map<pair<string, string >, size_t > nodepair_lid;

    vector<CellChat> CC_ve;

    // Deconvolution complex
    map<string, string > gene_complex;
    map<string, set<string > > complex_genes;

    // 

    GeneSet( string infile )
	{}
    GeneSet( )
	{}

    void readingeneset( string infile );
    void readinlinks( string infile );
    void readincellchat( string infile );
    void readincellchat3( string infile );
    void readindeconvolute( string infile );
    void readinoverlappedgenes( string infile );

    void deconvolutegetset();
    void AssignCellchat2link( Bank &bk );
    void AssignCellchat2link2( Bank &bk );

    void readinsigovlgenes( string infile, double pthr, int type );
    void readinsharedgokegg( string infile, set<string > &proi_go );
    void readincellchat2( string infile );

};


void Bank::readinconsgeneset( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
     while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );

        string ct = ps[0];
        string marker = ps[1];

        cons_geneset[ct].insert(marker);
    }
    inf.close();
}

void Bank::readinmodulegeneset( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );

        string gn = ps[0];
        string module = ps[5];

        module_geneset[module].insert(gn);
    }
    inf.close();
}

set<string > Bank::getgeneset(string item, string fea )
{
    if ( fea == "C" )
    {
        if ( cons_geneset.find( item ) == cons_geneset.end() )
        {
            cout<<"error cannot find item "<<item<<" in cons geneset "<<endl; exit(1);
        } 
        return cons_geneset[item];
    } else if ( fea == "M" )
    {
        if ( module_geneset.find( item ) == module_geneset.end() )
        {
            cout<<"error cannot find item "<<item<<" in module geneset "<<endl; exit(1);
        } 
        return module_geneset[item];
    } else{
        cout<<"error cannot identify fea ()"<<fea<<"()"<<endl; 
        cout<<"1\n";
        cout<<"2\n";
        exit(1);
    }
}

void GeneSet::readingeneset( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );

        string item = ps[0];
        string ct = ps[2];
        string fea = ps[3];
        items_celltype[item] = ct;
        items_feature[item] = fea;
    }
    inf.close();
}

void GeneSet::readinlinks( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string a = ps[0];
        string b = ps[1];
        double cor = atof(ps[2].c_str());
        double p = atof(ps[3].c_str() );
        Links lk(a, b, cor, p );
        link_ve.push_back( lk );
        nodepair_lid[make_pair(a,b)] = link_ve.size()-1;
    }
    inf.close();
}

void GeneSet::readinoverlappedgenes( string infile )
{
    map<pair<string, string >, set<string > > link_geneset;
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string nodea = ps[0];
        string nodeb = ps[1];
        string genes = ps[5];
        if ( genes == "NONE" )
            continue;
        set<string > geneset;
        vector<string > gene_ps = parse_string( genes, ',');
        
        for ( size_t i = 0; i < gene_ps.size(); ++i )
        {
            geneset.insert( gene_ps[i] );
        }
        link_geneset[make_pair(nodea, nodeb)] = geneset;
    }
    inf.close();

    for (size_t i = 0; i < link_ve.size(); ++i )
    {
        string nodea = link_ve[i].nodea;
        string nodeb = link_ve[i].nodeb;
        if ( link_geneset.find(make_pair(nodea, nodeb) ) != link_geneset.end() )
        {
            link_ve[i].overlapped_genes = link_geneset[make_pair(nodea, nodeb)];
        }
    }
}



void GeneSet::readincellchat( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string id = ps[0];
        CellChat cc(ps[0], ps[2], ps[3], ps[4], ps[5], ps[12], ps[20], ps[21], ps[17], ps[22], ps[23]);
        CC_ve.push_back(cc);
    }
    inf.close();
}

void GeneSet::readincellchat3( string infile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        CellChat cc(ps[0], ps[1]);
        CC_ve.push_back(cc);
    }
    inf.close();
}

void GeneSet::readindeconvolute(string infile )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
            break;
        vector<string > ps = parse_string(line);
        string gene = ps[0];
        string comp = ps[4];
        gene_complex[gene] = comp;
        complex_genes[comp].insert(gene);
    }
    inf.close();
}

set<string > deconvolute( string partner, string gene, map<string, set<string > >& complex_genes )
{
    set<string > genes;
    vector<string > ps = parse_string( partner, ':');
    if ( ps.size() < 2 )
    {
        cout<<"error unexpected partner "<<partner<<endl; exit(1);
    }
    if ( ps[0] == "simple" )
    {
        if ( gene == "" )
        {
            cout<<"error simple genes "<<endl; exit(1);

        }
        genes.insert( gene );
    } else if ( ps[0] == "complex")
    {
        string cp = ps[1];
        if ( ps.size() > 2 )
        {
            for ( size_t i = 2; i < ps.size(); ++i )
            {
                cp += (":"+ps[i]);
            }
        }
        if ( complex_genes.find( cp ) == complex_genes.end() )
        {
            cout<<"error cannot find complex_genes for "<<cp<<endl; exit(1);
        }
        genes = complex_genes[cp];
    } else{
        cout<<"error unexpected title for deconvolute "<<ps[0]<<endl; exit(1);
    }
    return genes;
}

void GeneSet::deconvolutegetset()
{
    for ( size_t i = 0; i < CC_ve.size(); ++i )
    {
        CC_ve[i].geneseta = deconvolute(CC_ve[i].partnera, CC_ve[i].genea, complex_genes );
        CC_ve[i].genesetb = deconvolute(CC_ve[i].partnerb, CC_ve[i].geneb, complex_genes );
        
    }
}

string obtaincelltype( string node, string fea )
{
    if ( fea == "M" )
    {
        vector<string > ps = parse_string( node, '_' );
        if ( (int)ps.size() <= 2 )
        {
            cout<<"error node name "<<node<<endl; exit(1);
        }
        string ct = "";
        for ( size_t i = 1; i < ps.size()-1; ++i )
        {
            
            if ( i == 1 )
            {
                ct = ps[i];
                
            } else{
                ct += ('_'+ps[i]);
            }
        } 
       
        return ct;
    } else if ( fea == "C" )
    {
        return node;
    } else
    {
        cout<<"unknown fea "<<fea<<endl; exit(1);
    }

}

bool checkgenein( set<string > &genea, set<string > &geneb )
{
    bool fd = false;
    for ( set<string >::iterator ite = genea.begin(); ite != genea.end(); ++ite )
    {
        if ( geneb.find( *ite ) != geneb.end() )
        {
            fd = true;
            break;
        }
    }
    return fd;
}

void GeneSet::AssignCellchat2link( Bank &bk )
{
    // first set each link for cell type and broad cell types
    map<pair<string, string >, set<size_t > > ct_linkid;
    for ( size_t i = 0; i < link_ve.size(); ++i )
    {
        string nodea = link_ve[i].nodea;
        string nodeb = link_ve[i].nodeb;
        link_ve[i].celltypea = obtaincelltype( nodea, items_feature[nodea] );
        link_ve[i].celltypeb = obtaincelltype( nodeb, items_feature[nodeb] );
        ct_linkid[make_pair(link_ve[i].celltypea, link_ve[i].celltypeb)].insert( i );
    //    cout<<link_ve[i].celltypea<<" "<<link_ve[i].celltypeb<<" "<<i<<endl;
    }

    for ( size_t i = 0; i < CC_ve.size(); ++i )
    {
        set<string > cc_genes_a = CC_ve[i].geneseta;
        set<string > cc_genes_b = CC_ve[i].genesetb;

        pair<string, string > p1 = make_pair( CC_ve[i].celltypea, CC_ve[i].celltypeb );
        pair<string, string > p2 = make_pair( CC_ve[i].broad_celltypea, CC_ve[i].broad_celltypeb );
        pair<string, string > p3 = make_pair( CC_ve[i].celltypea, CC_ve[i].broad_celltypeb );
        pair<string, string > p4 = make_pair( CC_ve[i].broad_celltypea, CC_ve[i].celltypeb );
        set<size_t > f_linkid;
        if ( ct_linkid.find( p1 ) != ct_linkid.end() )
        {
            f_linkid.insert( ct_linkid[p1].begin(), ct_linkid[p1].end() );
        }
        if ( ct_linkid.find( p2 ) != ct_linkid.end() )
        {
            f_linkid.insert( ct_linkid[p2].begin(), ct_linkid[p2].end() );
        }
        if ( ct_linkid.find( p3 ) != ct_linkid.end() )
        {
            f_linkid.insert( ct_linkid[p3].begin(), ct_linkid[p3].end() );
        }
        if ( ct_linkid.find( p4 ) != ct_linkid.end() )
        {
            f_linkid.insert( ct_linkid[p4].begin(), ct_linkid[p4].end() );
        }
        for ( set<size_t >::iterator ite = f_linkid.begin(); ite != f_linkid.end(); ++ite )
        {
            string itema = link_ve[*ite].nodea;
            string itemb = link_ve[*ite].nodeb;
            string feaa = items_feature[itema];
            string feab = items_feature[itemb];
            set<string > nodea_genes = bk.getgeneset( itema, feaa );
            set<string > nodeb_genes = bk.getgeneset( itemb, feab );
            bool fda = checkgenein( cc_genes_a, nodea_genes );
            bool fdb = checkgenein( cc_genes_b, nodeb_genes );
            if ( fda && fdb )
            {
                link_ve[*ite].CC_ind.insert( i );
                link_ve[*ite].CC_cps.insert( CC_ve[i].id);
            }
            
        }


        pair<string, string > p5 = make_pair( CC_ve[i].celltypeb, CC_ve[i].celltypea );
        pair<string, string > p6 = make_pair( CC_ve[i].broad_celltypeb, CC_ve[i].broad_celltypea );
        pair<string, string > p7 = make_pair( CC_ve[i].celltypeb, CC_ve[i].broad_celltypea );
        pair<string, string > p8 = make_pair( CC_ve[i].broad_celltypeb, CC_ve[i].celltypea );
        set<size_t > r_linkid;
        if ( ct_linkid.find( p5 ) != ct_linkid.end() )
        {
            r_linkid.insert( ct_linkid[p5].begin(), ct_linkid[p5].end() );
        }
        if ( ct_linkid.find( p6 ) != ct_linkid.end() )
        {
            r_linkid.insert( ct_linkid[p6].begin(), ct_linkid[p6].end() );
        }
        if ( ct_linkid.find( p7 ) != ct_linkid.end() )
        {
            r_linkid.insert( ct_linkid[p7].begin(), ct_linkid[p7].end() );
        }
        if ( ct_linkid.find( p8 ) != ct_linkid.end() )
        {
            r_linkid.insert( ct_linkid[p8].begin(), ct_linkid[p8].end() );
        }

        for ( set<size_t >::iterator ite = r_linkid.begin(); ite != r_linkid.end(); ++ite )
        {
            string itema = link_ve[*ite].nodea;
            string itemb = link_ve[*ite].nodeb;
            string feaa = items_feature[itema];
            string feab = items_feature[itemb];
            set<string > nodea_genes = bk.getgeneset( itema, feaa );
            set<string > nodeb_genes = bk.getgeneset( itemb, feab );
            bool fda = checkgenein( cc_genes_b, nodea_genes );
            bool fdb = checkgenein( cc_genes_a, nodeb_genes );
            if ( fda && fdb )
            {
                link_ve[*ite].CC_ind.insert( i );
                link_ve[*ite].CC_cps.insert( CC_ve[i].id);
            }
            
        }
    }
}

void GeneSet::AssignCellchat2link2( Bank &bk )
{
    
    for ( size_t i = 0; i < CC_ve.size(); ++i )
    {
        string ga = CC_ve[i].genea;
        string gb = CC_ve[i].geneb;
        
        for ( size_t j = 0; j < link_ve.size(); ++j )
        {
            string nodea = link_ve[j].nodea;
            string nodeb = link_ve[j].nodeb;
            string feaa = items_feature[nodea];
            string feab = items_feature[nodeb];
            set<string > nodea_genes = bk.getgeneset( nodea, feaa );
            set<string > nodeb_genes = bk.getgeneset( nodeb, feab );

            if ( nodea_genes.find( ga ) != nodea_genes.end() && nodeb_genes.find(gb) != nodeb_genes.end() )
            {
                link_ve[j].CC_ind.insert(i);
            } else if ( nodea_genes.find( gb ) != nodea_genes.end() && nodeb_genes.find(ga) != nodeb_genes.end() )
            {
                link_ve[j].CC_ind.insert(i);
            }

        }
    }


    
}

void GeneSet::readinsigovlgenes( string infile, double pthr, int type ) // type==1, inct; type==2, crct
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
            break;
        vector<string > ps = parse_string(line, '\t');
        string na = ps[0];
        string nb = ps[1];
        double p = atof(ps[7].c_str());
        if ( p <= pthr )
        {
            if ( nodepair_lid.find(make_pair(na, nb ) ) == nodepair_lid.end() )
            {
                cout<<"error cannot find pair "<<na<<" "<<nb<<" in links"<<endl; exit(1);
            }
            size_t lid = nodepair_lid[make_pair(na, nb )];
            if ( type == 1 )
                link_ve[lid].sig_ovl_gene_inct = true;
            else if ( type == 2)
                link_ve[lid].sig_ovl_gene_crct = true;
            else
            {
                cout<<"unkown type "<<type<<endl; exit(1);
            }
        }
    }
    inf.close();

}

void readinpriorgo( string infile, int nthr, set<string > &pgo )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
            break;
        vector<string > ps = parse_string(line, '\t');
        int n = atoi( ps[0].c_str() );
        string got = ps[1];
        if ( n >= nthr )
            pgo.insert( got );
    }
    inf.close();
}

void GeneSet::readinsharedgokegg( string infile, set<string > &proi_go )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
            break;
        vector<string > ps = parse_string(line, '\t');
        string na = ps[0];
        string nb = ps[1];
        if ( ps[2] == "NONE")
            continue;
        vector<string > gps = parse_string( ps[2], '|' );
        bool fd = false;
        for ( size_t i= 0; i < gps.size(); ++i )
        {
            string g = gps[i];
            if ( proi_go.find( g ) != proi_go.end() )
            {
                fd = true;
            }
        }
        if ( fd )
        {
            if ( nodepair_lid.find(make_pair(na, nb ) ) == nodepair_lid.end() )
            {
                cout<<"error cannot find pair "<<na<<" "<<nb<<" in links"<<endl; exit(1);
            }
            size_t lid = nodepair_lid[make_pair(na, nb )];
            link_ve[lid].shared_gokegg = true;
        }
    }
    inf.close();
}

void GeneSet::readincellchat2( string infile )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
            break;
        vector<string > ps = parse_string(line, '\t');
        string na = ps[0];
        string nb = ps[1];
        if ( nodepair_lid.find(make_pair(na, nb ) ) == nodepair_lid.end() )
        {
            cout<<"error cannot find pair "<<na<<" "<<nb<<" in links"<<endl; exit(1);
        }
        size_t lid = nodepair_lid[make_pair(na, nb )];
        link_ve[lid].CC = true;
    }
    inf.close();
}

void output_link_sumfea( GeneSet &gs, string outfile )
{
    ofstream outf( outfile.data() );

    outf<<"nodeA\tnodeB\tOvlGeneIn\tOvlGeneCr\tShareGO\tCC"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        outf<<nodea<<"\t"<<nodeb;
        if ( gs.link_ve[i].sig_ovl_gene_inct )
        {
            outf<<"\t1";
        } else{
            outf<<"\t0";
        }
        if ( gs.link_ve[i].sig_ovl_gene_crct )
        {
            outf<<"\t1";
        } else{
            outf<<"\t0";
        }
        if ( gs.link_ve[i].shared_gokegg )
        {
            outf<<"\t1";
        } else{
            outf<<"\t0";
        }
        if ( gs.link_ve[i].CC )
        {
            outf<<"\t1";
        } else{
            outf<<"\t0";
        }
        outf<<endl;
    }
    outf.close();
}

void getcombinedgenes( GeneSet &gs, Bank &bk, string outfile )
{
    set<string > cb_gene;
    for ( map<string, string >::iterator ite = gs.items_feature.begin(); 
        ite != gs.items_feature.end(); ++ite )
    {
        string item = ite->first;
        string fea = ite->second;
        if ( fea == "C" )
        {
            if ( bk.cons_geneset.find( item ) == bk.cons_geneset.end() )
            {
                cout<<"error cannot find item "<<item<<" in cons geneset "<<endl; exit(1);
            } else{
                cb_gene.insert( bk.cons_geneset[item].begin(), bk.cons_geneset[item].end() );
            }

        } else if ( fea == "M" )
        {
            if ( bk.module_geneset.find( item ) == bk.module_geneset.end() )
            {
                cout<<"error cannot find item "<<item<<" in module geneset "<<endl; exit(1);
            } else{
                cb_gene.insert( bk.module_geneset[item].begin(), bk.module_geneset[item].end() );
            }
        }

    }

    ofstream outf( outfile.data() );
    for ( set<string >::iterator ite = cb_gene.begin(); ite!= cb_gene.end(); ++ite )
    {
        outf<<*ite<<endl;
    }
    outf.close();
}

void cal_degree( GeneSet &gs, string outfile )
{
    map<string, map<string, vector< double > > > ct_ct_cors; 
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        string celltypea = gs.items_celltype[nodea];
        string celltypeb = gs.items_celltype[nodeb];
        double cor = gs.link_ve[i].cor;
        ct_ct_cors[celltypea][celltypeb].push_back(cor);
        ct_ct_cors[celltypeb][celltypea].push_back(cor);
    }
    ofstream outf( outfile.data() );
    for ( map<string, map<string, vector< double > > >::iterator ite = ct_ct_cors.begin(); ite != ct_ct_cors.end(); ++ite )
    {
        string cta = ite->first;
        outf<<">>"<<cta<<endl;
        double s = 0;
        for (map<string, vector< double > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string ctb = si->first;
            double a = 0;
            for ( vector<double >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
            {
                double c = *ti;
                if ( c > a )
                    a = c;
            }
            outf<<ctb<<"\t"<<a<<endl;
            s += a;
        }
        outf<<"Total\t"<<s<<endl;
        

    }
    outf<<endl;
}

void intersect_geneset( set<string > &seta, set<string > &setb, set<string > &itset )
{
    for ( set<string >::iterator ite = seta.begin(); ite != seta.end(); ++ite )
    {
        if ( setb.find(*ite ) != setb.end() )
            itset.insert( *ite );
    }
}

/*
void get_ovlgene_for_links( GeneSet &gs, Bank &bk, set<string > &all_conserved_marker, set<string> &all_module_gene, 
    string outfile, string outfile2, string outfile3, string outfile4, string outfile5 )
{

    ofstream outf( outfile.data() );

    outf<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP\tIntGenes"<<endl;

    ofstream outf3( outfile3.data() );

    outf3<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP\tIntGenes"<<endl;

    ofstream outf4( outfile4.data() );

    outf4<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP"<<endl;

    ofstream outf5( outfile5.data() );
    outf5<<"OvlGene\tOccur\tCelltypes"<<endl;

    set<string > totalintgene;
    map<string, set<string > > ovlgene_celltype;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
        }
        if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
        }
        string cta = gs.items_celltype[nodea];
        string ctb = gs.items_celltype[nodeb];
        string fea = gs.items_feature[nodea];
        string feb = gs.items_feature[nodeb];
        if ( cta != ctb )
        {
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > intset;
            intersect_geneset( genesetA, genesetB, intset );

            for ( set<string >::iterator si = intset.begin(); si != intset.end(); ++si )
            {
                ovlgene_celltype[*si].insert(cta);
                ovlgene_celltype[*si].insert(ctb);
            }
            totalintgene.insert( intset.begin(), intset.end() );
            int ca = (int)genesetA.size();
            int cb = (int)genesetB.size();
            int ci = (int)intset.size();

            if ( fea == feb )
            {
                if ( fea == "C" )
                {
                    int tt = (int)all_conserved_marker.size();
                    if ( ca >= cb )
                    {
                        double p = hypergeometrictest(ca, cb, tt, ci);
                        double fc = ( ci*1.0 / cb ) / ( ca*1.0 / tt );
                        outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                        if ( p < 0.05 )
                            outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    } else
                    {
                        double p = hypergeometrictest(cb, ca, tt, ci);
                        double fc = ( ci*1.0 / ca ) / ( cb*1.0 / tt );
                        outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                        if ( p < 0.05 )
                            outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    }
                } else if ( fea == "M" )
                {
                    int tt = (int)all_module_gene.size();
                    if ( ca >= cb )
                    {
                        double p = hypergeometrictest(ca, cb, tt, ci);
                        double fc = ( ci*1.0 / cb ) / ( ca*1.0 / tt );
                        outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                        if ( p < 0.05 )
                            outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    } else
                    {
                        double p = hypergeometrictest(cb, ca, tt, ci);
                        double fc = ( ci*1.0 / ca ) / ( cb*1.0 / tt );
                        outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                        if ( p < 0.05 )
                            outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    }
                }
            } else 
            {
                if ( fea == "C" && feb == "M" )
                {
                    int tt = (int)all_module_gene.size();
                    set<string > success;
                    intersect_geneset( genesetA, all_module_gene, success );
                    int cs = (int)success.size();
                    double p = 1;
                    double fc = 1;
                    if ( cs > 0 )
                    {
                        p = hypergeometrictest(cs, cb, tt, ci);
                        fc = ( ci*1.0 / cb ) / ( cs*1.0 / tt );

                    }
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                } else if ( fea == "M" && feb == "C" )
                {
                    int tt = (int)all_module_gene.size();
                    set<string > success;
                    intersect_geneset( genesetB, all_module_gene, success );
                    int cs = (int)success.size();
                    double p = 1;
                    double fc = 1;
                    if ( cs > 0 )
                    {
                        p = hypergeometrictest(cs, ca, tt, ci);
                        fc = ( ci*1.0 / ca ) / ( cs*1.0 / tt );

                    }
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                }
            }
             
            if ( intset.empty() )
            {
                outf<<"\tNONE"<<endl;
                
            } else
            {
                
                for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                {
                    if ( ite == intset.begin() )
                    {
                        outf<<"\t"<<*ite;
                    } else{
                        outf<<","<<*ite;
                    }
                }
                outf<<endl;
            }

        } else
        {
        
        }


    }
    outf.close();
    outf3.close();
    outf4.close();

    ofstream outf2(outfile2.data() );
    for ( set<string >::iterator ite = totalintgene.begin(); ite != totalintgene.end(); ++ite )
    {
        outf2<<*ite<<endl;
    }
    outf2.close();

    for ( map<string, set<string > >::iterator ite = ovlgene_celltype.begin(); ite != ovlgene_celltype.end(); ++ite )
    {
        string gn = ite->first;
        outf5<<gn<<"\t"<<(int)ite->second.size();
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            if ( si == ite->second.begin() )
                outf5<<"\t";
            else
                outf5<<",";
            outf5<<*si;
        }
        outf5<<endl;
    }
    outf5.close();
} */

void get_ovlgene_for_links_mc( GeneSet &gs, Bank &bk, set<string > &all_conserved_marker, set<string> &all_module_gene, 
    string outfile, string outfile2, string outfile3, string outfile4, string outfile5 )
{

    int sample_times = 1000;
    int cm_size = (int)all_conserved_marker.size();
    int mg_size = (int)all_module_gene.size();
    vector<string > acm_ve;
    for ( set<string >::iterator ite = all_conserved_marker.begin(); ite != all_conserved_marker.end(); ++ite )
    {
        acm_ve.push_back(*ite);
    }
    vector<string > amg_ve;
    for ( set<string >::iterator ite = all_module_gene.begin(); ite != all_module_gene.end(); ++ite )
    {
        amg_ve.push_back(*ite);
    }

    ofstream outf( outfile.data() );

    outf<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP\tIntGenes"<<endl;

    ofstream outf3( outfile3.data() );

    outf3<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP\tIntGenes"<<endl;

    ofstream outf4( outfile4.data() );

    outf4<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tfc\tOverlapP"<<endl;

    ofstream outf5( outfile5.data() );
    outf5<<"OvlGene\tOccur\tCelltypes"<<endl;

    set<string > totalintgene;
    map<string, set<string > > ovlgene_celltype;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        if ( i % 10 == 0 )
            cout<<i<<endl;
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
        }
        if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
        }
        string cta = gs.items_celltype[nodea];
        string ctb = gs.items_celltype[nodeb];
        string fea = gs.items_feature[nodea];
        string feb = gs.items_feature[nodeb];
        if ( cta != ctb )
        {
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > intset;
            intersect_geneset( genesetA, genesetB, intset );

            for ( set<string >::iterator si = intset.begin(); si != intset.end(); ++si )
            {
                ovlgene_celltype[*si].insert(cta);
                ovlgene_celltype[*si].insert(ctb);
            }
            totalintgene.insert( intset.begin(), intset.end() );
            int ca = (int)genesetA.size();
            int cb = (int)genesetB.size();
            int ci = (int)intset.size();

            // Monte-caro sampling 
            vector<int > sampled_intset_num;
            if ( fea == "C" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_intset;
                    intersect_geneset( s_genesetA, s_genesetB, s_intset );
                    sampled_intset_num.push_back((int)s_intset.size() );
                }
            } else if ( fea == "C" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_intset;
                    intersect_geneset( s_genesetA, s_genesetB, s_intset );
                    sampled_intset_num.push_back((int)s_intset.size() );
                }
            } else if ( fea == "M" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_intset;
                    intersect_geneset( s_genesetA, s_genesetB, s_intset );
                    sampled_intset_num.push_back((int)s_intset.size() );
                }
            } else if ( fea == "M" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_intset;
                    intersect_geneset( s_genesetA, s_genesetB, s_intset );
                    sampled_intset_num.push_back((int)s_intset.size() );
                }
            }
            double p = 0;
            double fc = 0;
            cal_mc_p_larger( (int)intset.size(), sampled_intset_num, p, fc );
            outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<fc<<"\t"<<p;
            if ( p < 0.05 )
                outf4<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<fc<<"\t"<<p<<endl;

            
             
            if ( intset.empty() )
            {
                outf<<"\tNONE"<<endl;
                
            } else
            {
                
                for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                {
                    if ( ite == intset.begin() )
                    {
                        outf<<"\t"<<*ite;
                    } else{
                        outf<<","<<*ite;
                    }
                }
                outf<<endl;
            }

        } else
        {
        /*  
            if ( fea != feb )
            {
                set<string > genesetA = bk.getgeneset(nodea, gs.items_feature[nodea] );
                set<string > genesetB = bk.getgeneset(nodeb, gs.items_feature[nodeb] );
                set<string > intset;
                intersect_geneset( genesetA, genesetB, intset );
                int cta = (int)genesetA.size();
                int ctb = (int)genesetB.size();
                int cti = (int)intset.size();
                double p = hypergeometrictest(cta, ctb, tt, cti);
                double oc = overlap_coefficient(cta, ctb, cti);
                outf3<<nodea<<"\t"<<nodeb<<"\t"<<cta<<"\t"<<ctb<<"\t"<<cti<<"\t"<<tt<<"\t"<<oc<<"\t"<<p;
                if ( intset.empty() )
                {
                    outf3<<"\tNONE"<<endl;
                } else
                {
                    if ( p < 0.01 )
                        outf4<<nodea<<"\t"<<nodeb<<"\t"<<cta<<"\t"<<ctb<<"\t"<<cti<<"\t"<<tt<<"\t"<<oc<<"\t"<<p<<endl;
            
                    for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                    {
                        if ( ite == intset.begin() )
                        {
                            outf3<<"\t"<<*ite;
                        } else{
                            outf3<<","<<*ite;
                        }
                    }
                    outf3<<endl;
                }
            }
            */
        }


    }
    outf.close();
    outf3.close();
    outf4.close();

    ofstream outf2(outfile2.data() );
    for ( set<string >::iterator ite = totalintgene.begin(); ite != totalintgene.end(); ++ite )
    {
        outf2<<*ite<<endl;
    }
    outf2.close();

    for ( map<string, set<string > >::iterator ite = ovlgene_celltype.begin(); ite != ovlgene_celltype.end(); ++ite )
    {
        string gn = ite->first;
        outf5<<gn<<"\t"<<(int)ite->second.size();
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            if ( si == ite->second.begin() )
                outf5<<"\t";
            else
                outf5<<",";
            outf5<<*si;
        }
        outf5<<endl;
    }
    outf5.close();
}

// for module and conserved
/*
void get_ovlgene_for_links2( GeneSet &gs, Bank &bk, int tt, string outfile )
{
    ofstream outf( outfile.data() );

    outf<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tOC\tOverlapP\tIntGenes"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
        }
        if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
        }
        string fea = gs.items_feature[nodea];
        string feb = gs.items_feature[nodeb];
        if ( fea != feb )
        {
            set<string > genesetA = bk.getgeneset(nodea, gs.items_feature[nodea] );
            set<string > genesetB = bk.getgeneset(nodeb, gs.items_feature[nodeb] );
            set<string > intset;
            intersect_geneset( genesetA, genesetB, intset );
            int ca = (int)genesetA.size();
            int cb = (int)genesetB.size();
            int ci = (int)intset.size();
            double p = hypergeometrictest(ca, cb, tt, ci);
            double oc = overlap_coefficient(ca, cb, ci);
            outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<oc<<"\t"<<p;

            if ( intset.empty() )
            {
                outf<<"\tNONE"<<endl;
            } else
            {
                for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                {
                    if ( ite == intset.begin() )
                    {
                        outf<<"\t"<<*ite;
                    } else{
                        outf<<","<<*ite;
                    }
                }
                outf<<endl;
            }

        }
    }
    outf.close();
}
*/

// for module across celltypes
/*
void get_ovlgene_for_links3( GeneSet &gs, Bank &bk, int tt, string outfile, string outfile2, string outfile3 )
{
    ofstream outf( outfile.data() );
    ofstream outf3( outfile3.data() );

    set<string > totalintgene;
    outf<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tOC\tOverlapP\tIntGenes"<<endl;
    outf3<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tIntCount\tGeneBankCount\tOC\tOverlapP\tIntGenes"<<endl;
    
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
        }
        if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
        }
        string cta = gs.items_celltype[nodea];
        string ctb = gs.items_celltype[nodeb];
        if ( cta != ctb )
        {
            set<string > genesetA = bk.getgeneset(nodea, gs.items_feature[nodea] );
            set<string > genesetB = bk.getgeneset(nodeb, gs.items_feature[nodeb] );
            set<string > intset;
            intersect_geneset( genesetA, genesetB, intset );
            totalintgene.insert( intset.begin(), intset.end() );
            int ca = (int)genesetA.size();
            int cb = (int)genesetB.size();
            int ci = (int)intset.size();
            double oc = overlap_coefficient( ca, cb, ci );
            double p = hypergeometrictest(ca, cb, tt, ci);
            outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<oc<<"\t"<<p;
            
            if ( intset.empty() )
            {
                outf<<"\tNONE"<<endl;
            } else
            {
                for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                {
                    if ( ite == intset.begin() )
                    {
                        outf<<"\t"<<*ite;
                    } else{
                        outf<<","<<*ite;
                    }
                }
                outf<<endl;
            }

            if ( p < 0.001 )
            {
                outf3<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<ci<<"\t"<<tt<<"\t"<<oc<<"\t"<<p;
                for ( set<string >::iterator ite = intset.begin(); ite != intset.end(); ++ite )
                {
                    if ( ite == intset.begin() )
                    {
                        outf3<<"\t"<<*ite;
                    } else{
                        outf3<<","<<*ite;
                    }
                }
                outf3<<endl;
            }

        }
    }
    outf.close();
    outf3.close();

    ofstream outf2(outfile2.data() );
    for ( set<string >::iterator ite = totalintgene.begin(); ite != totalintgene.end(); ++ite )
    {
        outf2<<*ite<<endl;
    }
    outf2.close();

}
*/

void get_lr_pair_for_gs( set<string > &gs, map<string, set<string > > &lrs, map<string, set<string > > &rls, set<string > &paired_gs )
{
    for ( set<string >::iterator ite = gs.begin(); ite != gs.end(); ++ite )
    {
        string g = *ite;
        if ( lrs.find(g) != lrs.end() )
        {
            paired_gs.insert(lrs[g].begin(), lrs[g].end() );
        }
        if ( rls.find(g) != rls.end() )
        {
            paired_gs.insert( rls[g].begin(), rls[g].end() );
        }
    }
}

void get_receptor_for_gs( set<string > &gs, map<string, set<string > > &lrs, set<string > &paired_gs )
{
    for ( set<string >::iterator ite = gs.begin(); ite != gs.end(); ++ite )
    {
        string g = *ite;
        if ( lrs.find(g) != lrs.end() )
        {
            paired_gs.insert(lrs[g].begin(), lrs[g].end() );
        }
        
    }
}

void get_ligand_for_gs( set<string > &gs, map<string, set<string > > &rls, set<string > &paired_gs )
{
    for ( set<string >::iterator ite = gs.begin(); ite != gs.end(); ++ite )
    {
        string g = *ite;
        if ( rls.find(g) != rls.end() )
        {
            paired_gs.insert(rls[g].begin(), rls[g].end() );
        }
        
    }
}

void get_lr_pair_for_pgs( set<string > &gs1, set<string > &gs2, map<string, set<string > > &lrs, map<string, set<string > > &rls,
    set<string> &in )
{
    
    for ( set<string >::iterator ite = gs1.begin(); ite != gs1.end(); ++ite )
    {
        string g = *ite;
        if ( lrs.find(g) != lrs.end() )
        {
            for ( set<string >::iterator si = lrs[g].begin(); si != lrs[g].end(); ++si )
            {
                if ( gs2.find(*si) != gs2.end() )
                {
                    string lrs_in = "f_"+g+"_"+*si;
                    in.insert(lrs_in);
                }
            }
        }
        if ( rls.find(g) != rls.end() )
        {
            for ( set<string >::iterator si = rls[g].begin(); si != rls[g].end(); ++si )
            {
                if ( gs2.find(*si) != gs2.end() )
                {
                    
                    string rls_in = "r_"+g+"_"+*si;
                    in.insert(rls_in);
                }
            }
        }
    }

}

void get_receptor_for_pgs( set<string > &gs1, set<string > &gs2, map<string, set<string > > &lrs, set<string > &in )
{
    for ( set<string >::iterator ite = gs1.begin(); ite != gs1.end(); ++ite )
    {
        string g = *ite;
        if ( lrs.find(g) != lrs.end() )
        {
            for ( set<string >::iterator si = lrs[g].begin(); si != lrs[g].end(); ++si )
            {
                if ( gs2.find(*si) != gs2.end() )
                {
                    string lrs_in = "f_"+g+"_"+*si;
                    in.insert(lrs_in);
                }
            }
        }
    }
}

void get_receptor_for_pgs( string cta, set<string > &gs1, string ctb, set<string > &gs2, 
    map<pair<string, string>, set<pair<string, string > > > &lrs, set<string > &in )
{
    for ( set<string >::iterator ite = gs1.begin(); ite != gs1.end(); ++ite )
    {
        string ga = *ite;
        pair<string, string > cga = make_pair(cta, ga);
        if ( lrs.find(cga) != lrs.end() )
        {
            for ( set<string >::iterator si = gs2.begin(); si != gs2.end(); ++si )
            {
                string gb = *si;
                pair<string, string > cgb = make_pair(ctb, gb);
                if ( lrs[cga].find(cgb) != lrs[cga].end() )
                {
                    string lrs_in = "f_"+ga+"_"+gb;
                    in.insert(lrs_in);
                }
            }
        }
    }
}


/*
void get_siglrp_for_links( GeneSet &gs, Bank &bk, 
    set<string > &all_conserved_marker, set<string> &all_module_gene, 
    map<string, set<string > > &lrs, map<string, set<string > > &rls, 
    string outfile1, string outfile2 )
{
    ofstream outf( outfile1.data() );
    ofstream outf2( outfile2.data() );
    outf<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tTotalSuc\tSampledSuc\tTotal\tfc\tOverlapP\tInLRP"<<endl;
    outf2<<"nodeA\tnodeB\tGeneCountA\tGeneGountB\tTotalSuc\tSampledSuc\tTotal\tfc\tOverlapP"<<endl;
    set<string > totalintgene;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
        }
        if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
        {
            cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
        }
        string cta = gs.items_celltype[nodea];
        string ctb = gs.items_celltype[nodeb];
        string fea = gs.items_feature[nodea];
        string feb = gs.items_feature[nodeb];
        if ( cta != ctb )
        {
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > gsA_pair;
            get_lr_pair_for_gs( genesetA, lrs, rls, gsA_pair );
            set<string > gsB_pair;
            get_lr_pair_for_gs( genesetB, lrs, rls, gsB_pair );
            set<string > success_gsA;
            set<string > success_gsB;
            if ( fea == "C" )
                intersect_geneset( all_conserved_marker, gsB_pair, success_gsA );
            else if ( fea == "M" )
                intersect_geneset( all_module_gene, gsB_pair, success_gsA );
            if ( feb == "C" )
                intersect_geneset( all_conserved_marker, gsA_pair, success_gsB );
            else if ( feb == "M")
                intersect_geneset( all_module_gene, gsA_pair, success_gsB );

            set<string > gsb_in_pair;
            set<string > gsa_in_pair;
            intersect_geneset( genesetB, gsA_pair, gsb_in_pair );
            intersect_geneset( genesetA, gsB_pair, gsa_in_pair );
            int ca = (int)genesetA.size();
            int cb = (int)genesetB.size();
            int csa = (int)success_gsA.size();
            int csb = (int)success_gsB.size();
            int ca_p = (int)gsa_in_pair.size();
            int cb_p = (int)gsb_in_pair.size();
            double p = 1;
            double fc = 1;

            

            if ( fea == feb )
            {
                int tt = 0;
                if ( fea == "C" )
                {
                    tt = (int)all_conserved_marker.size();
                    
                } else if ( fea == "M" )
                    tt = (int)all_module_gene.size();
                else
                {
                    cout<<"error cannot find right type "<<fea<<endl; exit(1);
                }

                if ( ca >= cb )
                {
                    p = hypergeometrictest(csb, cb, tt, cb_p );
                    if ( csb > 0 )
                        fc = ( cb_p*1.0 / cb ) / ( csb*1.0 / tt );
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csb<<"\t"<<cb_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf2<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csb<<"\t"<<cb_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    
                } else 
                {
                    p = hypergeometrictest(csa, ca, tt, ca_p );
                    if ( csa > 0 )
                        fc = ( ca_p*1.0 / ca ) / ( csa*1.0 / tt );
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csa<<"\t"<<ca_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf2<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csa<<"\t"<<ca_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    
                }
            } else 
            {   
                int tt = (int)all_module_gene.size();
                if ( fea == "C" && feb =="M")
                {
                    
                    p = hypergeometrictest(csb, cb, tt, cb_p );
                    if ( csb > 0 )
                        fc = ( cb_p*1.0 / cb ) / ( csb*1.0 / tt );
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csb<<"\t"<<cb_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf2<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csb<<"\t"<<cb_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    
                } else if ( fea == "M" && feb == "C" )
                {
                    p = hypergeometrictest(csa, ca, tt, ca_p );
                    if ( csa > 0 )
                        fc = ( ca_p*1.0 / ca ) / ( csa*1.0 / tt );
                    outf<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csa<<"\t"<<ca_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p;
                    if ( p < 0.05 )
                        outf2<<nodea<<"\t"<<nodeb<<"\t"<<ca<<"\t"<<cb<<"\t"<<csa<<"\t"<<ca_p<<"\t"<<tt<<"\t"<<fc<<"\t"<<p<<endl;
                    
                }
            }

            set<string > ins;
            
            get_lr_pair_for_pgs( genesetA, genesetB, lrs, rls, ins );
            if ( !ins.empty() )
            {
                for (set<string>::iterator ite = ins.begin(); ite != ins.end(); ++ite )
                {
                    if ( ite == ins.begin() )
                    {
                        outf<<"\t"<<*ite;
                    } else
                    {
                        outf<<"|"<<*ite;
                    }
                }
                outf<<endl;
            } else {
                outf<<"\tNONE"<<endl;
            }
        }
    }

    outf.close();
    outf2.close();
}
*/

void readinodnodes( string infile, vector<string > &nodes )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf,line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );

        string n = ps[0];
        nodes.push_back(n);
    }
    inf.close();
}
/*
void getlrt_matrix( GeneSet &gs, Bank &bk, vector<string > &nodes,
    set<string > &all_conserved_marker, set<string> &all_module_gene, 
    map<string, set<string > > &lrs, map<string, set<string > > &rls, 
    string outfile1, string outfile2, string outfile3 )
{
    map<string, map<string, double > > nodel_noder_p;
    map<string, map<string, set<string > > > nodel_noder_lrp;
    ofstream outf3( outfile3.data() );
    outf3<<"nodea\tnodeb\tCountA\tCountB\tCLa\tCRb\tSucc\tTotal\tP"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];

            if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
            }
            if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
            }
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            string fea = gs.items_feature[nodea];
            string feb = gs.items_feature[nodeb];
            if ( cta == ctb )
                continue;
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > gsA_receptor;
            get_receptor_for_gs( genesetA, lrs, gsA_receptor );
            set<string > success_gsB;
            if ( feb == "C" )
                intersect_geneset( all_conserved_marker, gsA_receptor, success_gsB );
            else if ( feb == "M")
                intersect_geneset( all_module_gene, gsA_receptor, success_gsB );
            set<string > gsB_ligand;
            get_ligand_for_gs( genesetB, rls, gsB_ligand );
            set<string > success_gsA;
            if ( fea == "C")
                intersect_geneset( all_conserved_marker, gsB_ligand, success_gsA );
            else if ( fea == "M")
                intersect_geneset( all_module_gene, gsB_ligand, success_gsA );


            set<string > gsb_in_pair;
            intersect_geneset( genesetB, gsA_receptor, gsb_in_pair );
            set<string > gsa_in_pair;
            intersect_geneset( genesetA, gsB_ligand, gsa_in_pair );

            outf3<<nodea<<"\t"<<nodeb<<"\t"<<(int)genesetA.size()<<"\t"<<(int)genesetB.size()<<"\t"<<(int)gsb_in_pair.size()<<"\t"<<(int)gsa_in_pair.size();
            if ( fea == feb )
            {
                int tt = 0;
                if ( feb == "C" )
                {
                    tt = (int)all_conserved_marker.size();
                    
                } else if ( feb == "M" )
                    tt = (int)all_module_gene.size();

                if ( cta >= ctb )
                {
                    int cb = (int)genesetB.size();
                    int csb = (int)success_gsB.size();
                    int cb_p = (int)gsb_in_pair.size();

                    double p = hypergeometrictest(csb, cb, tt, cb_p );
                    outf3<<"\t"<<csb<<"\t"<<tt<<"\t"<<p<<endl;
                    nodel_noder_p[nodea][nodeb] = p;
                } else 
                {
                    int ca = (int)genesetA.size();
                    int csa = (int)success_gsA.size();
                    int ca_p = (int)gsa_in_pair.size();
                    double p = hypergeometrictest(csa, ca, tt, ca_p );
                    outf3<<"\t"<<csa<<"\t"<<tt<<"\t"<<p<<endl;
                    nodel_noder_p[nodea][nodeb] = p;
                }
            } else{
            
                if ( fea == "C" && feb == "M" )
                {
                    int tt = (int)all_module_gene.size();
                    int cb = (int)genesetB.size();
                    int csb = (int)success_gsB.size();
                    int cb_p = (int)gsb_in_pair.size();
                    
                    double p = hypergeometrictest(csb, cb, tt, cb_p );
                    outf3<<"\t"<<csb<<"\t"<<tt<<"\t"<<p<<endl;
                    nodel_noder_p[nodea][nodeb] = p;
                } else if ( fea == "M" && feb == "C" )
                { 
                    int tt = (int)all_module_gene.size();
                    int ca = (int)genesetA.size();
                    int csa = (int)success_gsA.size();
                    int ca_p = (int)gsa_in_pair.size();
                    
                    double p = hypergeometrictest(csa, ca, tt, ca_p );
                    outf3<<"\t"<<csa<<"\t"<<tt<<"\t"<<p<<endl;
                    nodel_noder_p[nodea][nodeb] = p;
                }
            }
   

            set<string > ins;
            
            get_receptor_for_pgs( genesetA, genesetB, lrs, ins );
            if ( !ins.empty() )
            {
                nodel_noder_lrp[nodea][nodeb] = ins;
            }
        }
    }

    ofstream outf1( outfile1.data() );
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        outf1<<"\t"<<nodes[i];

    }
    outf1<<endl;
    ofstream outf2( outfile2.data() );
    outf2<<"nodea\tnodeb\tlrp"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        outf1<<nodea;
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            if ( cta == ctb )
            {
                outf1<<"\tNA";
            } else{
                double p = 1;
                if ( nodel_noder_p.find( nodea ) != nodel_noder_p.end() )
                {
                    if ( nodel_noder_p[nodea].find(nodeb) != nodel_noder_p[nodea].end() )
                    {
                        p = nodel_noder_p[nodea][nodeb];
                    }
                }
                double tp = 0;
                if ( p < 0.05 )
                    tp = (-1)*(log10(p+0.000001));
                outf1<<"\t"<<tp;

                set<string > lrps = nodel_noder_lrp[nodea][nodeb];
                if ( !lrps.empty() )
                {
                    outf2<<nodea<<"\t"<<nodeb;
                    for ( set<string >::iterator ite = lrps.begin(); ite != lrps.end(); ++ite )
                    {
                        if ( ite == lrps.begin() )
                        {
                            outf2<<"\t"<<*ite;
                        } else{
                            outf2<<"|"<<*ite;
                        }
                    }
                    outf2<<endl;
                }

            }
            
        }
        outf1<<endl;
    }
    outf1.close();
    outf2.close();
}
*/


void getlrt_matrix2( GeneSet &gs, Bank &bk, vector<string > &nodes,
    set<string > &all_conserved_marker, set<string> &all_module_gene, 
    map<string, set<string > > &lrs, map<string, set<string > > &rls, 
    string outfile1, string outfile2, string outfile3 )
{
    int sample_times = 1000;
    int cm_size = (int)all_conserved_marker.size();
    int mg_size = (int)all_module_gene.size();
    vector<string > acm_ve;
    for ( set<string >::iterator ite = all_conserved_marker.begin(); ite != all_conserved_marker.end(); ++ite )
    {
        acm_ve.push_back(*ite);
    }
    vector<string > amg_ve;
    for ( set<string >::iterator ite = all_module_gene.begin(); ite != all_module_gene.end(); ++ite )
    {
        amg_ve.push_back(*ite);
    }
    

    map<string, map<string, double > > nodel_noder_p;
    map<string, map<string, set<string > > > nodel_noder_lrp;
    ofstream outf3( outfile3.data() );
    outf3<<"nodea\tnodeb\tCountA\tCountB\tLRP_n\tLRP_fc\tP"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            cout<<i<<" "<<j<<endl;

            if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
            }
            if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
            }
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            string fea = gs.items_feature[nodea];
            string feb = gs.items_feature[nodeb];
            if ( cta == ctb )
                continue;
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > in_lrp;
            get_receptor_for_pgs( genesetA, genesetB, lrs, in_lrp );
            nodel_noder_lrp[nodea][nodeb] = in_lrp;

            // Monte-caro sampling 
            vector<int > sampled_lrp_num;
            if ( fea == "C" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( s_genesetA, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "C" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( s_genesetA, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "M" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( s_genesetA, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "M" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( s_genesetA, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            }
            double p = 0;
            double fc = 0;
            cal_mc_p_larger( (int)in_lrp.size(), sampled_lrp_num, p, fc );

            nodel_noder_p[nodea][nodeb] = p;

            outf3<<nodea<<"\t"<<nodeb<<"\t"<<(int)genesetA.size()<<"\t"<<(int)genesetB.size()<<"\t"<<(int)in_lrp.size()<<"\t"<<fc<<"\t"<<p<<endl;
        } 
    }

    ofstream outf1( outfile1.data() );
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        outf1<<"\t"<<nodes[i];

    }
    outf1<<endl;
    ofstream outf2( outfile2.data() );
    outf2<<"nodea\tnodeb\tlrp"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        outf1<<nodea;
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            if ( cta == ctb )
            {
                outf1<<"\tNA";
            } else{
                double p = 1;
                if ( nodel_noder_p.find( nodea ) != nodel_noder_p.end() )
                {
                    if ( nodel_noder_p[nodea].find(nodeb) != nodel_noder_p[nodea].end() )
                    {
                        p = nodel_noder_p[nodea][nodeb];
                    }
                }
                double tp = 0;
                if ( p < 0.05 )
                    tp = (-1)*(log10(p+0.000001));
                outf1<<"\t"<<tp;

                set<string > lrps = nodel_noder_lrp[nodea][nodeb];
                if ( !lrps.empty() )
                {
                    outf2<<nodea<<"\t"<<nodeb;
                    for ( set<string >::iterator ite = lrps.begin(); ite != lrps.end(); ++ite )
                    {
                        if ( ite == lrps.begin() )
                        {
                            outf2<<"\t"<<*ite;
                        } else{
                            outf2<<"|"<<*ite;
                        }
                    }
                    outf2<<endl;
                }

            }
            
        }
        outf1<<endl;
    }
    outf1.close();
    outf2.close();

}

void getlrt_matrix2( GeneSet &gs, Bank &bk, vector<string > &nodes,
    set<string > &all_conserved_marker, set<string> &all_module_gene, 
    map<pair<string, string>, set<pair<string, string > > > &lrs, 
    map<pair<string, string>, set<pair<string, string > > > &rls, 
    string outfile1, string outfile2, string outfile3 )
{
    int sample_times = 1000;
    int cm_size = (int)all_conserved_marker.size();
    int mg_size = (int)all_module_gene.size();
    vector<string > acm_ve;
    for ( set<string >::iterator ite = all_conserved_marker.begin(); ite != all_conserved_marker.end(); ++ite )
    {
        acm_ve.push_back(*ite);
    }
    vector<string > amg_ve;
    for ( set<string >::iterator ite = all_module_gene.begin(); ite != all_module_gene.end(); ++ite )
    {
        amg_ve.push_back(*ite);
    }
    

    map<string, map<string, double > > nodel_noder_p;
    map<string, map<string, set<string > > > nodel_noder_lrp;
    ofstream outf3( outfile3.data() );
    outf3<<"nodea\tnodeb\tCountA\tCountB\tLRP_n\tLRP_fc\tP"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            
            cout<<i<<" "<<j<<endl;

            if ( gs.items_celltype.find( nodea ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodea "<<nodea<<" in items_celltype "<<endl; exit(1);
            }
            if ( gs.items_celltype.find( nodeb ) == gs.items_celltype.end() )
            {
                cout<<"error cannot find nodeb "<<nodeb<<" in items_celltype "<<endl; exit(1);
            }
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            string fea = gs.items_feature[nodea];
            string feb = gs.items_feature[nodeb];
            if ( cta == ctb )
                continue;
            set<string > genesetA = bk.getgeneset(nodea, fea );
            set<string > genesetB = bk.getgeneset(nodeb, feb );
            set<string > in_lrp;
            get_receptor_for_pgs( cta, genesetA, ctb, genesetB, lrs, in_lrp );
            nodel_noder_lrp[nodea][nodeb] = in_lrp;

            // Monte-caro sampling 
            vector<int > sampled_lrp_num;
            if ( fea == "C" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( cta, s_genesetA, ctb, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "C" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( cm_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( acm_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( cta, s_genesetA, ctb, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "M" && feb == "C" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( cm_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( acm_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( cta, s_genesetA, ctb, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            } else if ( fea == "M" && feb == "M" )
            {
                vector<vector<int > > shaffled_index_A;
                shaffle_vector( mg_size, (int)genesetA.size(), sample_times, shaffled_index_A );
                vector<vector<int > > shaffled_index_B;
                shaffle_vector( mg_size, (int)genesetB.size(), sample_times, shaffled_index_B );
                for ( size_t si = 0; si < shaffled_index_A.size(); ++si )
                {
                    set<string > s_genesetA;
                    get_sampled_geneset( amg_ve, shaffled_index_A[si], s_genesetA );
                    set<string > s_genesetB;
                    get_sampled_geneset( amg_ve, shaffled_index_B[si], s_genesetB );
                    set<string > s_lrp;
                    get_receptor_for_pgs( cta, s_genesetA, ctb, s_genesetB, lrs, s_lrp );
                    sampled_lrp_num.push_back((int)s_lrp.size() );
                }
            }
            double p = 0;
            double fc = 0;
            cal_mc_p_larger( (int)in_lrp.size(), sampled_lrp_num, p, fc );

            nodel_noder_p[nodea][nodeb] = p;

            outf3<<nodea<<"\t"<<nodeb<<"\t"<<(int)genesetA.size()<<"\t"<<(int)genesetB.size()<<"\t"<<(int)in_lrp.size()<<"\t"<<fc<<"\t"<<p<<endl;
        } 
    }

    ofstream outf1( outfile1.data() );
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        outf1<<"\t"<<nodes[i];

    }
    outf1<<endl;
    ofstream outf2( outfile2.data() );
    outf2<<"nodea\tnodeb\tlrp"<<endl;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        outf1<<nodea;
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            if ( cta == ctb )
            {
                outf1<<"\tNA";
            } else{
                double p = 1;
                if ( nodel_noder_p.find( nodea ) != nodel_noder_p.end() )
                {
                    if ( nodel_noder_p[nodea].find(nodeb) != nodel_noder_p[nodea].end() )
                    {
                        p = nodel_noder_p[nodea][nodeb];
                    }
                }
                double tp = 0;
                if ( p < 0.05 )
                    tp = (-1)*(log10(p+0.000001));
                outf1<<"\t"<<tp;

                set<string > lrps = nodel_noder_lrp[nodea][nodeb];
                if ( !lrps.empty() )
                {
                    outf2<<nodea<<"\t"<<nodeb;
                    for ( set<string >::iterator ite = lrps.begin(); ite != lrps.end(); ++ite )
                    {
                        if ( ite == lrps.begin() )
                        {
                            outf2<<"\t"<<*ite;
                        } else{
                            outf2<<"|"<<*ite;
                        }
                    }
                    outf2<<endl;
                }

            }
            
        }
        outf1<<endl;
    }
    outf1.close();
    outf2.close();

}


void output_link_with_CC( GeneSet &gs, string outfile, string outfile2 )
{
    ofstream outf( outfile.data() );
    outf<<"source\ttarget\tCCsize\tCC_cp"<<endl;
    map<string, set<string > > cp_pair;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        
        if ( gs.link_ve[i].CC_cps.empty() )
        {
        //    outf<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t0\tNONE"<<endl;
        } else
        {
            outf<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t"<<gs.link_ve[i].CC_cps.size();
            string ps = gs.link_ve[i].nodea+":"+gs.link_ve[i].nodeb;
            for ( set<string >::iterator ite = gs.link_ve[i].CC_cps.begin();
                ite != gs.link_ve[i].CC_cps.end(); ++ite )
            {
                cp_pair[*ite].insert(ps);
                if ( ite == gs.link_ve[i].CC_cps.begin() )
                {
                    outf<<"\t"<<*ite;
                } else{
                    outf<<","<<*ite;
                }
            }
            outf<<endl;
        }
    }
    outf.close();

    ofstream outf2( outfile2.data() );
    for ( map<string, set<string > >::iterator ite = cp_pair.begin(); ite != cp_pair.end(); ++ite )
    {
        string cp = ite->first;
        outf2<<cp;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            if ( si == ite->second.begin() )
            {
                outf2<<"\t"<<*si;
            } else{
                outf2<<","<<*si;
            }
        }
        outf2<<endl;
    }
    outf2.close();
}

void output_link_with_CC2( GeneSet &gs, string outfile )
{
    ofstream outf( outfile.data() );
    outf<<"source\ttarget\tCCsize"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        
        if ( gs.link_ve[i].CC_ind.empty() )
        {
        //    outf<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t0"<<endl;
        } else
        {
            outf<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t"<<gs.link_ve[i].CC_ind.size();
            
            outf<<endl;
        }
    }
    outf.close();
}

void get_go_kegg( string infile, map<string, set<string > > &term_genes )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string term = ps[1];
        string genestr = ps[5];

        vector<string > gene_ps = parse_string( genestr, ',');
        for ( size_t i = 0; i < gene_ps.size(); ++i )
        {
            term_genes[term].insert( gene_ps[i] );
        }
    }
    inf.close();
}

void ana_overlapped_go_kegg( Bank &bk, GeneSet &gs, map<string, set<string > > &term_genes,
    string outfile1, string outfile2 )
{
    // mode1: require overlapped genes in term
    ofstream outf1(outfile1.data() );
    outf1<<"nodaA\tnodeB\tovlgeneCount\tovlgeneinTerm"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        if ( gs.link_ve[i].overlapped_genes.empty() )
            continue;
        
        map<string, int > term_gc;
        for ( map<string, set<string > >::iterator ti = term_genes.begin(); 
            ti != term_genes.end(); ++ti )
        {
            set<string > ovlgenes;
            intersect_geneset( gs.link_ve[i].overlapped_genes, ti->second, ovlgenes );
            if ( !ovlgenes.empty() )
            {
                term_gc[ti->first] = (int)ovlgenes.size();
            }
        }
        if ( term_gc.empty() )
            continue;

        outf1<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t"<<(int)gs.link_ve[i].overlapped_genes.size();
        for ( map<string, int >::iterator ti = term_gc.begin(); ti != term_gc.end(); ++ti )
        {
            string term = ti->first;
            int c = ti->second;
            if ( ti == term_gc.begin() )
            {
                outf1<<"\t"<<term<<":"<<c;
            } else{
                outf1<<"|"<<term<<":"<<c;
            }
            
        }
        outf1<<endl;
    }
    outf1.close();

    // mode2 require both 20% genes in shared term 
    ofstream outf2( outfile2.data() );
    outf2<<"nodaA\tnodeB\tovlgeneCount\tsharedTerm"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        set<string > nodea_gene = bk.getgeneset(nodea, gs.items_feature[nodea]);
        set<string > nodeb_gene = bk.getgeneset(nodeb, gs.items_feature[nodeb]);
        map<string, pair<int, int > > term_ca_cb;
        for ( map<string, set<string > >::iterator ti = term_genes.begin(); 
            ti != term_genes.end(); ++ti )
        {
            string term = ti->first;
            set<string > ovl_genea;
            intersect_geneset( nodea_gene, ti->second, ovl_genea );
            set<string > ovl_geneb;
            intersect_geneset( nodeb_gene, ti->second, ovl_geneb );

            if ( ovl_genea.empty() || ovl_geneb.empty() )
                continue;

            int cna = (int)nodea_gene.size();
            int coa = (int)ovl_genea.size();
            int cnb = (int)nodeb_gene.size();
            int cob = (int)ovl_geneb.size();
            int ctg = (int)ti->second.size();
            double fna = (coa*1.0)/cna;
            double fta = (coa*1.0)/ctg;
            double fnb = (cob*1.0)/cnb;
            double ftb = (cob*1.0)/ctg;
            if ( (fna > 0.1 || fta > 0.1 ) && (fnb > 0.1 || ftb > 0.1) )
            {
                term_ca_cb[term] = make_pair(coa, cob );
            }
        }

        if ( term_ca_cb.empty() )
        {
            continue;
        }

        outf2<<gs.link_ve[i].nodea<<"\t"<<gs.link_ve[i].nodeb<<"\t"<<(int)gs.link_ve[i].overlapped_genes.size();
        for ( map<string, pair<int, int > >::iterator ti = term_ca_cb.begin(); ti != term_ca_cb.end(); ++ti )
        {
            string term = ti->first;
            
            if ( ti == term_ca_cb.begin() )
            {
                outf2<<"\t"<<term<<":"<<ti->second.first<<":"<<ti->second.second;
            } else{
                outf2<<"|"<<term<<":"<<ti->second.first<<":"<<ti->second.second;
            }
            
        }
        outf2<<endl;

    }
    outf2.close();

}

void get_go_kegg2( string infile, map<string, set<string > > &term_gokegg )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string term = infile.substr(0, infile.find(".GOKEGG"));
    cout<<term<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string t = ps[1];
        term_gokegg[term].insert( t );
    }
    inf.close();
}

void get_go_kegg2_b( string infile, map<string, set<string > > &term_gokegg )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        get_go_kegg2( line, term_gokegg );
    }
    inf.close();
}


void ana_overlapped_go_kegg2( Bank &bk, GeneSet &gs, map<string, set<string > > &term_gokegg,
    string outfile, string outfile2 )
{
    ofstream outf1(outfile.data() );
    outf1<<"nodaA\tnodeB\tovlTerm"<<endl;
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        if ( term_gokegg.find( nodea ) == term_gokegg.end() || term_gokegg.find( nodeb ) == term_gokegg.end() )
        {
            outf1<<nodea<<"\t"<<nodeb<<"\tNONE"<<endl;
            continue;
        }
        set<string > gka = term_gokegg[nodea];
        set<string > gkb = term_gokegg[nodeb];
        set<string > ovlgk;
        intersect_geneset(gka, gkb, ovlgk );
        if ( ovlgk.empty() )
        {
            outf1<<nodea<<"\t"<<nodeb<<"\tNONE"<<endl;
            continue;
        }
        outf1<<nodea<<"\t"<<nodeb;
        for ( set<string >::iterator ite = ovlgk.begin(); ite != ovlgk.end(); ++ite )
        {
            if ( ite == ovlgk.begin() )
            {
                outf1<<"\t"<<*ite;
            } else{
                outf1<<"|"<<*ite;
            }
        }
        outf1<<endl;
    }
    outf1.close();

    map<string, set<string > > gk_ct;
    map<string, set<string > > gk_term;
    for ( map<string, set<string > >::iterator ite = term_gokegg.begin(); ite != term_gokegg.end(); ++ite )
    {
        string term = ite->first;
        string ct = gs.items_celltype[term];
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string g = *si;
            gk_ct[g].insert( ct );
            gk_term[g].insert( term );
        }
    }
    multimap<int, string, greater<int> > order_gk;
    for ( map<string, set<string > >::iterator ite = gk_ct.begin(); ite != gk_ct.end(); ++ite )
    {
        string gk = ite->first;
        int n = ite->second.size();
        order_gk.insert( make_pair(n, gk) );
    }
    ofstream outf2( outfile2.data() );
    for ( multimap<int, string, greater<int> >::iterator ite = order_gk.begin(); ite != order_gk.end(); ++ite )
    {
        outf2<<ite->first<<"\t"<<ite->second;
        for ( set<string >::iterator si = gk_term[ite->second].begin(); si != gk_term[ite->second].end(); ++si )
        {
            if ( si == gk_term[ite->second].begin() )
            {
                outf2<<"\t"<<*si;
            } else{
                outf2<<","<<*si;
            }
            
        }
        outf2<<endl;
    }
    outf2.close();

}

void readinkeygoterms( string infile, map<string, set<string > > &term_genes  )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string goid = ps[0];
        string gene_str = ps[10];
        vector<string > pgs = parse_string( gene_str, ',');
        set<string > genes;
        for ( size_t i = 0; i < pgs.size(); ++i )
        {
            genes.insert( pgs[i] );
        }
        term_genes[goid] = genes;
    }
    inf.close();
}

void ana_overlapped_keygo( Bank &bk, GeneSet &gs, set<string > &term_genes,
    string outfile )
{
    int k = 0;
    for ( set<string >::iterator ite = term_genes.begin(); ite != term_genes.end(); ++ite )
    {
        cout<<*ite<<endl;
        k += 1;
        if ( k == 10 )
            break;
    }
    ofstream outf( outfile.data() );
    outf<<"nodeA\tnodeB\tGeneCountA\tGeneCountB\tIntermA\tIntermB\tRatioA\tRatioB\tIntermA_gene\tIntermB_gene\n";
    for ( size_t i = 0; i < gs.link_ve.size(); ++i )
    {
        string nodea = gs.link_ve[i].nodea;
        string nodeb = gs.link_ve[i].nodeb;
        string celltypea = gs.items_celltype[nodea];
        string celltypeb = gs.items_celltype[nodeb];
        
        if ( celltypea == celltypeb )
            continue;
        set<string > genes_a = bk.getgeneset(nodea, gs.items_feature[nodea]);
        set<string > genes_b = bk.getgeneset(nodeb, gs.items_feature[nodeb]);
        set<string > genes_a_interm;
        intersect_geneset( genes_a, term_genes, genes_a_interm );
        set<string > genes_b_interm;
        intersect_geneset( genes_b, term_genes, genes_b_interm );

        outf<<nodea<<"\t"<<nodeb<<"\t"<<(int)genes_a.size()<<"\t"<<(int)genes_b.size();
        outf<<"\t"<<(int)genes_a_interm.size()<<"\t"<<(int)genes_b_interm.size();
        double ratioa = (1.0*(int)genes_a_interm.size())/(int)genes_a.size();
        double ratiob = (1.0*(int)genes_b_interm.size())/(int)genes_b.size();
        outf<<"\t"<<ratioa<<"\t"<<ratiob;
        for ( set<string >::iterator si = genes_a_interm.begin(); si != genes_a_interm.end(); ++si )
        {
            if ( si == genes_a_interm.begin() )
                outf<<"\t"<<*si;
            else
                outf<<","<<*si;
        }
       
        for ( set<string >::iterator si = genes_b_interm.begin(); si != genes_b_interm.end(); ++si )
        {
            if ( si == genes_b_interm.begin() )
                outf<<"\t"<<*si;
            else
                outf<<","<<*si;
        }
        outf<<endl;

    }
    outf.close();
}

void readinkeygokegg( string infile, set<string > &keygk )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        keygk.insert(line);
    }
    inf.close();
}

void ana_shared_gokegg( set<string > &keygk, string infile, string outfile )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line );
    ofstream outf( outfile.data() );
    outf<<"nodaA\tnodeB\tCount\tKey\n";
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string nodea = ps[0];
        string nodeb = ps[1];
        string term = ps[2];
        if ( term == "NONE")
            continue;
        vector<string > tps = parse_string( term, '|');
        int ts = (int)tps.size();
        bool fd = false;
        for ( size_t i = 0; i < tps.size(); ++i )
        {
            if ( keygk.find( tps[i]) != keygk.end() )
                fd = true;
        }
        outf<<nodea<<"\t"<<nodeb<<"\t"<<ts;
        if ( fd )
            outf<<"\t1"<<endl;
        else
            outf<<"\t0"<<endl;
    }
    outf.close();
}

void readinligandtarget( string infile, string gm, map<string, map<string, set<string > > > &l_t_map )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line );
    
    while(!inf.eof())
	{
        getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, ',' );
        string ligand = ps[0];
        string target = ps[1];
        l_t_map[gm][ligand].insert( target );

    }
    inf.close();
}

void readinligandtarget_b( string infile, map<string, map<string, set<string > > > &l_t_map )
{
    ifstream inf( infile.data() );
    if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    
    while(!inf.eof())
	{
        getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
        string filename = ps[0];
        string gs = ps[1];
        readinligandtarget( filename, gs, l_t_map );
    }
    inf.close();
}

pair<bool, bool> check_ligand_in_lrp( set<pair<string, string > > &sender_receiv, 
    GeneSet &gs, vector<string > &nodes, string nodea, string nodeb )
{
    string cta = gs.items_celltype[nodea];
    string ctb = gs.items_celltype[nodeb];
    bool suc = false;
    bool suc_in_a = false;
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        string node = nodes[i];
        string ct = gs.items_celltype[node];
        if ( ct == ctb )
            continue;
        if ( sender_receiv.find( make_pair(ct, ctb ) ) != sender_receiv.end() )
        {
            suc = true;
            break;
        }
    }
    if ( sender_receiv.find( make_pair(cta, ctb ) ) != sender_receiv.end() )
    {
        suc_in_a = true;
    }

    return make_pair(suc, suc_in_a );
}

class DOT_g
{
public:
    map<string, set<string > > ligand_target;
    set<string > ligands;
    set<string > targets;

    map<string, set<string > > ligand_node;
    map<string, set<string > > target_node;

    DOT_g( )
	{}
    
    void generate_dot_g( Bank &bk, GeneSet &gs, string outfile );

};

void DOT_g::generate_dot_g( Bank &bk, GeneSet &gs, string outfile )
{
    ofstream outf( outfile.data() );
    outf<<"digraph G {"<<endl;
    outf<<"\trankdir=LR; ranksep=.75; nodesep=.1; size=\"7.5,7.5\""<<endl;
    set<string > ligand_node_comb;
    map<string, set<string > > node_ligand;
    for ( map<string, set<string > >::iterator ite = ligand_node.begin(); ite != ligand_node.end(); ++ite )
    {
        ligand_node_comb.insert( ite->second.begin(), ite->second.end() );
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            node_ligand[*si].insert( ite->first );
        }
    }
    set<string > target_node_comb;
    map<string, set<string > > node_target;
    for ( map<string, set<string > >::iterator ite = target_node.begin(); ite != target_node.end(); ++ite )
    {
        target_node_comb.insert( ite->second.begin(), ite->second.end() );
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            node_target[*si].insert( ite->first );
        }
    }
    set<string > ligand_celltype;
    map<string, set<string > > celltype_ligandnode;
    for ( set<string >::iterator ite = ligand_node_comb.begin(); ite != ligand_node_comb.end(); ++ite )
    {
        string ct = gs.items_celltype[*ite];
        ligand_celltype.insert( ct );
        celltype_ligandnode[ct].insert(*ite);
    }
    set<string > target_celltype;
    map<string, set<string > > celltype_targetnode;
    for ( set<string >::iterator ite = target_node_comb.begin(); ite != target_node_comb.end(); ++ite )
    {
        string ct = gs.items_celltype[*ite];
        target_celltype.insert( ct );
        celltype_targetnode[ct].insert(*ite);
    }

    outf<<"\t{"<<endl;
    outf<<"\t\tnode [shape=plaintext];"<<endl;
    outf<<"\t\t\"Broad_cell_type\" -> \"Gene_set\" -> \"Ligand\" -> \"Target\" -> \" Gene_set\" -> \" Broad_cell_type\" [arrowhead=none];"<<endl;
    outf<<"\t}"<<endl;

    outf<<"\tnode [shape=plaintext];"<<endl;
    outf<<"\t{ rank = same; \"Broad_cell_type\"; ";
    for ( set<string >::iterator ite = ligand_celltype.begin(); ite != ligand_celltype.end(); ++ite )
    {
        outf<<"\":"<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\t{ rank = same; \"Gene_set\"; ";
    for ( set<string >::iterator ite = ligand_node_comb.begin(); ite != ligand_node_comb.end(); ++ite )
    {
        outf<<"\""<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\tnode [shape=plaintext];"<<endl;
    outf<<"\t{ rank = same; \"Ligand\"; ";
    for ( set<string >::iterator ite = ligands.begin(); ite != ligands.end(); ++ite )
    {
        outf<<"\""<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\t{ rank = same; \"Target\"; ";
    for ( set<string >::iterator ite = targets.begin(); ite != targets.end(); ++ite )
    {
        outf<<"\" "<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\tnode [shape=plaintext];"<<endl;
   
    outf<<"\t{ rank = same; \" Gene_set\"; ";
    for ( set<string >::iterator ite = target_node_comb.begin(); ite != target_node_comb.end(); ++ite )
    {
        outf<<"\" "<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\t{ rank = same; \" Broad_cell_type\"; ";
    for ( set<string >::iterator ite = target_celltype.begin(); ite != target_celltype.end(); ++ite )
    {
        outf<<"\" :"<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;

    for ( map<string, set<string > >::iterator ite = celltype_ligandnode.begin(); ite != celltype_ligandnode.end(); ++ite )
    {
        string ct = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string node = *si;
            outf<<"\t\":"<<ct<<"\" -> \""<<node<<"\" [arrowhead=none];"<<endl;
        }
    }
    for ( map<string, set<string > >::iterator ite = node_ligand.begin(); ite != node_ligand.end(); ++ite )
    {
        string node = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string ligand = *si;
            outf<<"\t\""<<node<<"\" -> \""<<ligand<<"\" [arrowhead=none,color=cyan];"<<endl;
        }
    }
    for ( map<string, set<string > >::iterator ite = ligand_target.begin(); ite != ligand_target.end(); ++ite )
    {
        string l = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string t = *si;
            outf<<"\t\""<<l<<"\" -> \" "<<t<<"\" [color=red];"<<endl;
        }
    }
    for ( map<string, set<string > >::iterator ite = node_target.begin(); ite != node_target.end(); ++ite )
    {
        string node = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string t = *si;
            outf<<"\t\" "<<t<<"\" -> \" "<<node<<"\" [arrowhead=none,color=cyan];"<<endl;
        }
    }
    for ( map<string, set<string > >::iterator ite = celltype_targetnode.begin(); ite != celltype_targetnode.end(); ++ite )
    {
        string ct = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string node = *si;
            outf<<"\t\" "<<node<<"\" -> \" :"<<ct<<"\" [arrowhead=none];"<<endl;
        }
    }
    outf<<"}"<<endl;
}

bool check_DOT_g_overlap( DOT_g &g1, DOT_g &g2, int thr_n )  // check target overlap
{
    bool ovl = false;
    int n_o = 0;
    for ( set<string >::iterator ite = g2.targets.begin(); ite != g2.targets.end(); ++ite )
    {
        if ( g1.targets.find(*ite) != g1.targets.end() )
        {
            
            n_o += 1;
        }
    }
   
    if ( n_o > thr_n )
        ovl = true;
    return ovl;
}

void merge_two_DOT_g( DOT_g &g1, DOT_g &g2 )  // put g2 into g1, and update g1
{
    for ( map<string, set<string > >::iterator ite = g2.ligand_target.begin(); ite != g2.ligand_target.end(); ++ite )
    {
        string l = ite->first;
        if ( g1.ligand_target.find( l ) != g1.ligand_target.end() )
        {
            g1.ligand_target[l].insert( ite->second.begin(), ite->second.end() );
        } else{
            g1.ligand_target.insert( *ite );
        }
    }
    for ( map<string, set<string > >::iterator ite = g2.ligand_node.begin(); ite != g2.ligand_node.end(); ++ite )
    {
        string l = ite->first;
        if ( g1.ligand_node.find( l ) != g1.ligand_node.end() )
        {
            g1.ligand_node[l].insert( ite->second.begin(), ite->second.end() );
        } else{
            g1.ligand_node.insert( *ite );
        }
    }
    for ( map<string, set<string > >::iterator ite = g2.target_node.begin(); ite != g2.target_node.end(); ++ite )
    {
        string l = ite->first;
        if ( g1.target_node.find( l ) != g1.target_node.end() )
        {
            g1.target_node[l].insert( ite->second.begin(), ite->second.end() );
        } else{
            g1.target_node.insert( *ite );
        }
    }
    g1.ligands.insert( g2.ligands.begin(), g2.ligands.end() );
    g1.targets.insert( g2.targets.begin(), g2.targets.end() );

}

void merge_DOT_g( vector<DOT_g > &sorted_DOT_g, vector<DOT_g> &merged_DOT_g, int thr )
{
    set<size_t > out_ids;
    for ( size_t i = 0; i < sorted_DOT_g.size(); ++i )
    {
        if ( out_ids.find(i) != out_ids.end() )
            continue;
        DOT_g new_g;
        new_g = sorted_DOT_g[i];
        out_ids.insert(i);
        for ( size_t j = i+1; j < sorted_DOT_g.size(); ++j )
        {
            if ( out_ids.find(j ) != out_ids.end() )
                continue;
            if ( check_DOT_g_overlap(new_g, sorted_DOT_g[j], thr ) )
            {
                merge_two_DOT_g( new_g, sorted_DOT_g[j] );
                out_ids.insert( j );
            }
        }
        merged_DOT_g.push_back(new_g );
    }
}

void generate_dotfile( vector<DOT_g > &merged_DOT_g, string prefix, Bank &bk, GeneSet &gs )
{
    size_t thr = 10;
    for ( size_t i = 0; i < merged_DOT_g.size(); ++i )
    {
        string outfile = prefix+".G"+inttostr((int)i)+".dot";
        merged_DOT_g[i].generate_dot_g( bk, gs, outfile );
    }
}

void generate_dot_sender_receiv( Bank &bk, GeneSet &gs, map<string, map<string, int > > &sr, string outprefix )
{
    string outfile = outprefix + ".sender_receiv.dot";
    ofstream outf( outfile.data() );
    outf<<"digraph G {"<<endl;
    outf<<"\trankdir=LR; ranksep=.75; nodesep=.1; size=\"7.5,7.5\""<<endl;
   
    set<string > senders;
    set<string > receivers;
    for ( map<string, map<string, int > >::iterator ite = sr.begin(); ite != sr.end(); ++ite )
    {
        string s = ite->first;
        senders.insert( s );
        for ( map<string, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            receivers.insert( si->first );
        }
    }

    set<string > sender_celltype;
    map<string, set<string > > celltype_sender;
    for ( set<string >::iterator ite = senders.begin(); ite != senders.end(); ++ite )
    {
        string ct = gs.items_celltype[*ite];
        sender_celltype.insert( ct );
        celltype_sender[ct].insert(*ite);
    }
    set<string > receiv_celltype;
    map<string, set<string > > celltype_receiv;
    for ( set<string >::iterator ite = receivers.begin(); ite != receivers.end(); ++ite )
    {
        string ct = gs.items_celltype[*ite];
        receiv_celltype.insert( ct );
        celltype_receiv[ct].insert(*ite);
    }

    outf<<"\t{"<<endl;
    outf<<"\t\tnode [shape=plaintext];"<<endl;
    outf<<"\t\t\"Broad_cell_type\" -> \"Sender\" -> \"Receiver\" -> \" Broad_cell_type\" [arrowhead=none];"<<endl;
    outf<<"\t}"<<endl;

    outf<<"\tnode [shape=plaintext];"<<endl;
    outf<<"\t{ rank = same; \"Broad_cell_type\"; ";
    for ( set<string >::iterator ite = sender_celltype.begin(); ite != sender_celltype.end(); ++ite )
    {
        outf<<"\":"<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\t{ rank = same; \"Sender\"; ";
    for ( set<string >::iterator ite = senders.begin(); ite != senders.end(); ++ite )
    {
        outf<<"\""<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;

    outf<<"\t{ rank = same; \"Receiver\"; ";
    for ( set<string >::iterator ite = receivers.begin(); ite != receivers.end(); ++ite )
    {
        outf<<"\" "<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;
    outf<<"\t{ rank = same; \" Broad_cell_type\"; ";
    for ( set<string >::iterator ite = receiv_celltype.begin(); ite != receiv_celltype.end(); ++ite )
    {
        outf<<"\" :"<<*ite<<"\"; ";
    }
    outf<<"}"<<endl;

    for ( map<string, set<string > >::iterator ite = celltype_sender.begin(); ite != celltype_sender.end(); ++ite )
    {
        string ct = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string node = *si;
            outf<<"\t\":"<<ct<<"\" -> \""<<node<<"\" [arrowhead=none];"<<endl;
        }
    }

    for ( map<string, map<string, int > >::iterator ite = sr.begin(); ite != sr.end(); ++ite )
    {
        string s = ite->first;
        senders.insert( s );
        for ( map<string, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            outf<<"\t\""<<s<<"\" -> \" "<<si->first<<"\" [weight="<<si->second<<"];"<<endl;
        }
    }
    for ( map<string, set<string > >::iterator ite = celltype_receiv.begin(); ite != celltype_receiv.end(); ++ite )
    {
        string ct = ite->first;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string node = *si;
            outf<<"\t\" "<<node<<"\" -> \" :"<<ct<<"\" [arrowhead=none];"<<endl;
        }
    }
    outf<<"}"<<endl;

}

void ana_ligand_target( Bank &bk, GeneSet &gs, vector<string > &nodes,
    set<string > &all_conserved_marker, set<string> &all_module_gene, 
    map<string, map<string, set<string > > > &l_t_map,
    map<string, set<pair<string, string > > > &ligand_sender_receiver,
    string outfile1, string outfile2, string outfile3, string prefix )
{

    ofstream outf1(outfile1.data() );    // link_l_t_p
    outf1<<"nodea\tnodeb\tn_ligand\tn_genesa\tn_target\ta_genesb\tratio_a\tratio_b\tnet"<<endl;

    map<string, set<string > > ligand_target;
    map<string, set<string > > target_ligand;
    map<string, set<string > > ligand_node;
    map<string, set<string > > target_node;
    set<string > combined_all_genes;
    combined_all_genes = all_conserved_marker;
    map<string, set<string > > node_target;
    map<string, set<string > > node_ligand;
    
    combined_all_genes.insert( all_module_gene.begin(), all_module_gene.end() );
    int n_allgene = (int)combined_all_genes.size();

    map<string, set<string> > celltype_map;   // broad cell type : subtypes map
    for ( map<string, string >::iterator ite = gs.items_celltype.begin(); ite != gs.items_celltype.end(); ++ite )
    {
        celltype_map[ite->second].insert(ite->first);
    }

    map<string, map<string, int > > sender_receiv_count;

    for (size_t i = 0; i < nodes.size(); ++i )
    {
        string nodea = nodes[i];
        set<string > genes_a = bk.getgeneset(nodea, gs.items_feature[nodea]);

        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            string nodeb = nodes[j];
            if ( j == i )
                continue;

            string cta = gs.items_celltype[nodea];
            string ctb = gs.items_celltype[nodeb];
            if ( cta == ctb )
                continue;

            set<string > genes_b = bk.getgeneset(nodeb, gs.items_feature[nodeb]);

            if ( l_t_map.find( nodeb) == l_t_map.end() )
            {
             //   cout<<"warning: no record for nodeb "<<nodeb<<endl; 
                outf1<<nodea<<"\t"<<nodeb<<"\t"<<0<<"\t"<<(int)genes_a.size()<<"\t"<<0<<"\t"<<(int)genes_b.size()<<"\t"<<0<<"\t"<<0<<"\tNA"<<endl;
                continue;
            }

            set<string > suc_ligand;
            set<string > suc_ligand_in_a;
            set<string > suc_target_in_b;

            int n = 0;

            for ( map<string, set<string > >::iterator ite = l_t_map[nodeb].begin(); ite != l_t_map[nodeb].end(); ++ite )
            {
                string ligand = ite->first;
                if ( combined_all_genes.find( ligand ) == combined_all_genes.end() )
                    continue;

                // check lrp
                if ( ligand_sender_receiver.find(ligand) == ligand_sender_receiver.end())
                    continue;

                pair<bool, bool> ligand_suc = check_ligand_in_lrp( ligand_sender_receiver[ligand], gs, nodes, nodea, nodeb );
                bool ligand_suc_a = ligand_suc.first;
                if ( !ligand_suc_a )
                {
                    continue;
                }
                suc_ligand.insert( ligand );

                
                if ( genes_a.find( ligand ) != genes_a.end() )
                {
                    bool ligand_suc_in_a = ligand_suc.second;
                    if ( ligand_suc_in_a )
                    {
                        suc_ligand_in_a.insert( ligand );
                        ligand_node[ligand].insert( nodea );
                        node_ligand[nodea].insert(ligand);

                        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
                        {
                            if ( genes_b.find(*si ) != genes_b.end() )
                            {
                                target_node[*si].insert( nodeb );
                                ligand_target[ligand].insert( *si );
                                target_ligand[*si].insert( ligand );
                                node_target[nodeb].insert(*si);
                                suc_target_in_b.insert(*si);
                                n += 1;
                            }
                            
                        }
                    }
                }
            }

            if ( n > 0 )
                sender_receiv_count[nodea][nodeb] = n;

 /*           if ( !suc_ligand.empty() )
            {
// hypergeometric test
// r: number of success
// n: sample size
// N: Total object
// k: number of success events drawn from n samples 
                

                int n_suc = (int)suc_ligand.size();
                int n_genesa = (int)genes_a.size();
                int n_suc_in_a = (int)suc_ligand_in_a.size();
                double P = hypergeometrictest(n_suc, n_genesa, n_allgene, n_suc_in_a );
                outf1<<nodea<<"\t"<<nodeb<<"\t"<<n_suc<<"\t"<<n_genesa<<"\t"<<n_allgene<<"\t"<<n_suc_in_a<<"\t"<<P<<endl;
            } */
            
            if ( !suc_ligand_in_a.empty() )
            {
                int n_ligand = (int)suc_ligand_in_a.size();
                int n_target = (int)suc_target_in_b.size();
                
                double ratio_a = n_ligand*1.0/(int)genes_a.size();
                double ratio_b = n_target*1.0/(int)genes_b.size();
                outf1<<nodea<<"\t"<<nodeb<<"\t"<<n_ligand<<"\t"<<(int)genes_a.size()<<"\t"<<n_target<<"\t"<<(int)genes_b.size()<<"\t"<<ratio_a<<"\t"<<ratio_b;
                for ( set<string >::iterator ti = suc_ligand_in_a.begin(); ti != suc_ligand_in_a.end(); ++ti ) 
                {

                    string l = *ti;
                    if ( ti == suc_ligand_in_a.begin() )
                    {
                        outf1<<"\t"<<l;
                    } else{
                        outf1<<"|"<<l;
                    }
                    for ( set<string >::iterator gi = suc_target_in_b.begin(); gi != suc_target_in_b.end(); ++gi )
                    {

                        if ( gi == suc_target_in_b.begin() )
                        {
                            outf1<<":"<<*gi;
                        } else
                        {
                            outf1<<","<<*gi;
                        }
                    }
                    
                }
                outf1<<endl;
            } else
            {
                outf1<<nodea<<"\t"<<nodeb<<"\t"<<0<<"\t"<<(int)genes_a.size()<<"\t"<<0<<"\t"<<(int)genes_b.size()<<"\t"<<0<<"\t"<<0<<"\tNA"<<endl;
            }

        }
    }

    // generete sender receiver map
//    generate_dot_sender_receiv( bk, gs, sender_receiv_count, prefix );

    multimap<int, string, greater<int > > size_ligand;
    for ( map<string, set<string > >::iterator ite = ligand_target.begin(); ite != ligand_target.end(); ++ite )
    {
        string l = ite->first;
        int n = 0;
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string t = *si;
            set<string > cts;
            for ( set<string >::iterator ti = target_node[t].begin(); ti != target_node[t].end(); ++ti )
            {
                string ct = gs.items_celltype[*ti];
                cts.insert(ct);
            }
            int cn = (int)cts.size();
            n += cn;
        }
        size_ligand.insert( make_pair(n, l) );
    }

    vector<DOT_g > sorted_DOT_g;
    
    string outfile5 = prefix+".sorted_target.stat.txt";
    ofstream outf5(outfile5.data() );
    outf5<<"Ligand\tn_sender\tn_target\tnreceiver"<<endl;
    ofstream outf2(outfile2.data() );
    for ( multimap<int, string, greater<int > >::iterator ite = size_ligand.begin(); ite != size_ligand.end(); ++ite )
    {
        string l = ite->second;

        DOT_g g;
        g.ligands.insert( l);
        outf2<<">"<<l;
        
        set<string> sender_cts;
        for ( set<string >::iterator si = ligand_node[l].begin(); si != ligand_node[l].end(); ++si )
        {
            if ( si == ligand_node[l].begin() )
            {
                outf2<<"\t";
            } else
                outf2<<":";
            outf2<<*si;

            g.ligand_node[l].insert(*si);
            string ct = gs.items_celltype[*si];
            sender_cts.insert(ct);
        }
        outf2<<endl;
        outf5<<l<<"\t"<<(int)sender_cts.size()<<"\t"<<(int)ligand_target[l].size();
        set<string > receivers_cts;
        for ( set<string >::iterator si = ligand_target[l].begin(); si != ligand_target[l].end(); ++si )
        {
            string t = *si;
            g.targets.insert( t);
            g.ligand_target[l].insert(t);
            outf2<<t;
            
            for ( set<string >::iterator ti = target_node[t].begin(); ti != target_node[t].end(); ++ti )
            {
                if ( ti == target_node[t].begin() )
                {
                    outf2<<"\t";
                } else
                    outf2<<":";
                outf2<<*ti;
                g.target_node[t].insert(*ti);
                string ct = gs.items_celltype[*ti];
                receivers_cts.insert(ct);
            }
            outf2<<endl;

        }
        outf5<<"\t"<<(int)receivers_cts.size()<<endl;

        sorted_DOT_g.push_back( g );

    }
    outf2.close();

    ofstream outf3( outfile3.data() );
    for ( map<string, set<string > >::iterator ite = node_target.begin(); ite != node_target.end(); ++ite )
    {
        set<string > genes_a = bk.getgeneset(ite->first, gs.items_feature[ite->first]);
        int genen = (int)genes_a.size();
        int nt = (int)ite->second.size();
        double ratio = nt*1.0/genen;
        outf3<<ite->first<<"\t"<<nt<<"\t"<<genen<<"\t"<<ratio<<endl;
    }
    outf3.close();

    string outfile4 = prefix+".ligand_count.txt";
    ofstream outf4( outfile4.data() );
    outf4<<"Node\tLigand_count"<<endl;
    
    for ( map<string, set<string > >::iterator ite = node_ligand.begin(); ite != node_ligand.end(); ++ite )
    {
        outf4<<ite->first<<"\t"<<(int)ite->second.size()<<endl;
    }
    outf4.close();
/*
    // generate DOT graph for ligand
    vector<DOT_g> merged_DOT_g;
    int ovl_n_thr = 5;
    merge_DOT_g( sorted_DOT_g, merged_DOT_g, ovl_n_thr );

    cout<<"generate_dotfile"<<endl;
    generate_dotfile( merged_DOT_g, prefix, bk, gs );
*/
}

void exit_with_help()
{
	cerr <<"integrate gene program and analysis"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-C		cons_celltype_marker file input" <<endl;
	cerr <<"-M		module gene file input" <<endl;
	cerr <<"-N		Group Node file input"<<endl;
	cerr <<"-L		Group Link file input"<<endl;
    cerr <<"-v      Link overlapped gene file infput"<<endl;
    cerr <<"-G      GOfile input"<<endl;
    cerr <<"-K      KEGGfile input"<<endl;
    cerr <<"-g      key GO file input"<<endl;
    cerr <<"-O      output prefix"<<endl;
    cerr <<"-c      CellphoneDB_result file input OR simplified CC result for addi==1"<<endl;
    cerr <<"-d      CellphoneDB_deconvolute file input"<<endl;
    cerr <<"-f      runfunction type [combinegene | linkgeneoverlap | assigncellchat | gokegg | keygo | summarize | sharedgokegg | siglrp | ligandtarget | connectiondegree ]"<<endl;
    cerr <<"-a      additional modes [1,2 for linkgeneoverlap; 1,2 for gokegg; 1 for assigncellchat]"<<endl;
    cerr <<"-i      link with overlapgene in cell type input"<<endl;
    cerr <<"-r      link with overlapgene across cell type input"<<endl;
    cerr <<"-p      prioritized go pathway file input"<<endl;
    cerr <<"-n      conserved_marker_gene_union file input"<<endl;
    cerr <<"-m      module gene union file input"<<endl;
    cerr <<"-l      combined l-r-p file input"<<endl;
    cerr <<"-t      combined ligend target input"<<endl;
    cerr <<"/nIf additional modes 1 for gokegg, the GOfile input indicate the bunch of GOKEGGfile"<<endl;
    exit(1);
    
}

void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}


int main( int argc, char* argv[] )
{
    // test
 /*  func1();
    cout<<"sec"<<endl;
    func1();
    return 1; */

    if (argc == 1)
	{
		exit_with_help();
	}
	
    string inconsgenefile = "";
    string inmodulegenefile = "";
    string innodesfile = "";
    string inlinksfile = "";
    string ccresultfile = "";
    string ccdeconvfile = "";
    string func = "";
    string outprefix = "out";
    string gofile = "";
    string keggfile = "";
    string keygofile = "";
    string linkoverlappedgenefile = "";
    string addi = "";
    string inovlgenefile="";
    string crovlgenefile ="";
    string priorgofile="";
    string consmarkerUfile="";
    string modulegeneUfile="";
    string lrpfile= "";
    string ltfile = "";

	for(int i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			exit_with_help("Options must start with \'-\'.");

		if(argv[i][2] != '\0')
			exit_with_help("The option should be exactly one letter.");
		int option = argv[i][1];

		i++;

		switch(option)
		{
		
		case 'C':
			inconsgenefile = argv[i];
			break;
        case 'M':
			inmodulegenefile = argv[i];
			break;
        case 'N':
			innodesfile = argv[i];
			break;
        case 'L':
            inlinksfile = argv[i];
            break;
        case 'v':
            linkoverlappedgenefile = argv[i];
            break;
        case 'f':
            func = argv[i];
            break;
        case 'O':
            outprefix = argv[i];
            break;
        case 'G':
            gofile = argv[i];
            break;
        case 'g':
            keygofile = argv[i];
            break;
        case 'K':
            keggfile = argv[i];
            break;
        case 'c':
            ccresultfile = argv[i];
            break;
        case 'd':
            ccdeconvfile = argv[i];
            break;
        case 'a':
            addi = argv[i];
            break;
        case 'i':
            inovlgenefile = argv[i];
            break;
        case 'r':
            crovlgenefile = argv[i];
            break;
        case 'p':
            priorgofile = argv[i];
            break;
        case 'm':
            modulegeneUfile = argv[i];
            break;
        case 'n':
            consmarkerUfile = argv[i];
            break;
        case 'l':
            lrpfile = argv[i];
            break;
        case 't':
            ltfile = argv[i];
            break;
        
        default:
			exit_with_help();
        }
    }

    if ( func == "combinegene" )
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        string outfile = outprefix+".combinedgene.txt";
        getcombinedgenes( gs, bk, outfile );

        return 1;

    } else if ( func == "connectiondegree" )
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );


        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        string outfile = outprefix+".Conndegree.txt";
        cal_degree( gs, outfile);

    }
    else if ( func == "linkgeneoverlap" )
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        if ( addi == "" )
        {
            if ( consmarkerUfile == "")
            {
                exit_with_help("Error No consmarkerUfile indicated!");
            }
            if ( modulegeneUfile == "" )
            {
                exit_with_help("Error No modulegeneUfile indicated!");
            }
            string outfile = outprefix+".link_overlapped_gene.txt";
            string outfile2 = outprefix+".link_overlapped_gene.combined.txt";
            string outfile3 = outprefix+".link_CM_overlapped_gene.txt";
            string outfile4 = outprefix+".link_sig_overlappedgene.txt";
            string outfile5 = outprefix+".ovlgene_celltypes.txt";

            set<string > consmarker_gs;
            set<string > modulegene_gs;
            readingeneset( consmarkerUfile, consmarker_gs );
            readingeneset( modulegeneUfile, modulegene_gs );
            
            cout<<"linkgeneoverlap"<<endl;

            get_ovlgene_for_links_mc( gs, bk, consmarker_gs, modulegene_gs, outfile, outfile2, outfile3, outfile4, outfile5 );

        } else if ( addi == "1" )
        {
            string outfile = outprefix+".link_overlapped_gene.txt";
            int tt = 30000;
        //    get_ovlgene_for_links2( gs, bk, tt, outfile );
            cout<<"Not valid any more for addi == 1"<<endl;
        } else if ( addi == "2")
        {
            
            string outfile = outprefix+".link_overlapped_gene.txt";
            string outfile2 = outprefix+".link_overlapped_gene.combined.txt";
            string outfile3 = outprefix+".link_overlapped_gene.sig.txt";
            int tt = 30000;
        //    get_ovlgene_for_links3( gs, bk, tt, outfile, outfile2, outfile3 );
            cout<<"Not valid any more for addi == 2"<<endl;
        }

        return 1;
    } else if ( func == "assigncellchat" )
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }


        if ( ccresultfile == "")
        {
            exit_with_help("Error No ccresultfile indicated!");
        }
        
        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        if ( addi == "" )
        {
            if ( ccdeconvfile == "" )
            {
                exit_with_help("Error No ccdeconvfile indicated!");
            }
            gs.readincellchat( ccresultfile );
            gs.readindeconvolute( ccdeconvfile );
            
            gs.deconvolutegetset();

            gs.AssignCellchat2link( bk );

            string outfile = outprefix+ ".link_with_CC.txt";
            string outfile2 = outprefix+".CC.txt";
            output_link_with_CC( gs, outfile, outfile2 );
        } else if ( addi == "1" )
        {
            gs.readincellchat3( ccresultfile );
            gs.AssignCellchat2link2( bk );

            string outfile = outprefix+".link_with_CC.s.txt";
            output_link_with_CC2( gs, outfile );
        }

        


    } else if ( func == "gokegg" ) 
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }

        if ( addi == "" )
        {
            if ( linkoverlappedgenefile == "" )
            {
                exit_with_help("Error No linkoverlappedgenefile indicated!");
            }
            if ( gofile == "" )
            {
                exit_with_help("Error No gofile indicated!");
            }
            if ( keggfile == "" )
            {
                exit_with_help("Error No keggfile indicated!");
            }
        } else if ( addi == "1" )
        {
            if ( gofile == "" )
            {
                exit_with_help("Error No gofile indicated!");
            }
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        if ( addi == "" )
        {
            gs.readinoverlappedgenes( linkoverlappedgenefile );

            map<string, set<string > > term_genes;
            get_go_kegg( gofile, term_genes );
            get_go_kegg( keggfile, term_genes );

            string outfile1 = outprefix + ".link_with_gokegg_mode1.txt";
            string outfile2 = outprefix + ".link_with_gokegg_mode2.txt";
            ana_overlapped_go_kegg( bk, gs, term_genes, outfile1, outfile2 );
        } else if ( addi == "1")
        {
            map<string, set<string > > term_gokegg;
            get_go_kegg2_b( gofile, term_gokegg );

            string outfile = outprefix + ".link_with_ovlgokegg.txt";
            string outfile2 = outprefix+".gokegg_prior.txt";
            ana_overlapped_go_kegg2( bk, gs, term_gokegg, outfile, outfile2 );
        }
        
    } else if ( func == "keygo" ) 
    {

        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }
        if ( keygofile == "" )
        {
            exit_with_help("Error No keygofile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        map<string, set<string > > term_genes;
        readinkeygoterms( keygofile, term_genes );
        
        string anato_dev_id = "GO:0048856";
        string outfile1 = outprefix+".link_with_anat_dev.txt";
        ana_overlapped_keygo( bk, gs, term_genes[anato_dev_id], outfile1 );

        string resp_id = "GO:0050896";
        string outfile2 = outprefix+".link_with_resp_stim.txt";
        ana_overlapped_keygo( bk, gs, term_genes[resp_id], outfile2 );

    } else if ( func == "summarize" ) 
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }
        if ( inovlgenefile == "" )
        {
            exit_with_help("Error No inovlgenefile indicated!");
        }
        if ( crovlgenefile == "" )
        {
            exit_with_help("Error No crovlgenefile indicated!");
        }
        if ( gofile == "" )
        {
            exit_with_help("Error No gofile indicated!");
        }
        if ( ccresultfile == "" )
        {
            exit_with_help("Error No ccresultfile indicated!");
        }
        if ( priorgofile == "" )
        {
            exit_with_help("Error No priorgofile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        double pthr = 0.01;
        gs.readinsigovlgenes( inovlgenefile, pthr, 1 );
        gs.readinsigovlgenes( crovlgenefile, pthr, 2 );

        int nthr = 2;
        set<string > pgo;
        readinpriorgo( priorgofile, nthr, pgo );

        gs.readinsharedgokegg( gofile, pgo );

        gs.readincellchat2( ccresultfile );

        string outfile = outprefix+".link_with_summary.txt";
        output_link_sumfea( gs, outfile );
        
    } else if ( func == "sharedgokegg" ) 
    {
        if ( gofile == "" )
        {
            exit_with_help("Error No linkwithgofile indicated!");
        } 
        if ( keygofile == "" )
        {
            exit_with_help("Error No keygofile indicated!");
        }

        set<string > keygokegg;
        readinkeygokegg( keygofile, keygokegg );

        string outfile = outprefix+".link_with_sharedgokegg.txt";
        ana_shared_gokegg( keygokegg, gofile,  outfile );


    } else if ( func == "siglrp" )
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }
        if ( consmarkerUfile == "")
        {
            exit_with_help("Error No consmarkerUfile indicated!");
        }
        if ( modulegeneUfile == "" )
        {
            exit_with_help("Error No modulegeneUfile indicated!");
        }
        if ( lrpfile == "" )
        {
            exit_with_help("Error No lrpfile indicated!");
        }
        string outfile1 = outprefix+".link_LRP.txt";
        string outfile2 = outprefix+".link_sigLRP.txt";

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        set<string > consmarker_gs;
        set<string > modulegene_gs;
        readingeneset( consmarkerUfile, consmarker_gs );
        readingeneset( modulegeneUfile, modulegene_gs );

    /*    map<string, set<string > > lrs; 
        map<string, set<string > > rls;
        readin_lrp( lrpfile, lrs, rls ); */
        map<pair<string, string >, set<pair<string, string > > > lrs;
        map<pair<string, string >, set<pair<string, string > > > rls;
        readin_lrp( lrpfile, lrs, rls );

    //    get_siglrp_for_links( gs, bk, consmarker_gs, modulegene_gs, lrs, rls, outfile1, outfile2 );
        
        vector<string > odnodes;
        readinodnodes( innodesfile , odnodes );

        string outfile3 = outprefix+".L2R.P_matrix.txt";
        string outfile4 = outprefix+".L2R.list.txt";
        string outfile5 = outprefix+".link_L2R.txt";
        getlrt_matrix2( gs, bk, odnodes, consmarker_gs, modulegene_gs,  lrs, rls, outfile3, outfile4, outfile5 );

    } else if ( func == "ligandtarget" ) 
    {
        if ( inconsgenefile == "" )
        {
            exit_with_help("Error No consgenefile indicated!");

        }
        if ( innodesfile == "")
        {
            exit_with_help("Error No innodesfile indicated!");
        }
        if ( inmodulegenefile == "")
        {
            exit_with_help("Error No inmodulegenefile indicated!");
        }
        if ( inlinksfile == "" )
        {
            exit_with_help("Error No inlinksfile indicated!");
        }
        if ( consmarkerUfile == "")
        {
            exit_with_help("Error No consmarkerUfile indicated!");
        }
        if ( modulegeneUfile == "" )
        {
            exit_with_help("Error No modulegeneUfile indicated!");
        }
        if ( ltfile == "" )
        {
            exit_with_help("Error No ltfile indicated!");
        }
        if ( lrpfile == "" )
        {
            exit_with_help("Error No lrpfile indicated!");
        }

        Bank bk;
        bk.readinconsgeneset( inconsgenefile );

        bk.readinmodulegeneset( inmodulegenefile );

        GeneSet gs;
        gs.readingeneset( innodesfile );

        gs.readinlinks( inlinksfile );

        set<string > consmarker_gs;
        set<string > modulegene_gs;
        readingeneset( consmarkerUfile, consmarker_gs );
        readingeneset( modulegeneUfile, modulegene_gs );

        map<string, map<string, set<string > > > l_t_map;
        readinligandtarget_b( ltfile, l_t_map );

        map<string, set<pair<string, string > > > ligand_sender_receiver;
        readin_lrp2( lrpfile, ligand_sender_receiver );

        vector<string > odnodes;
        readinodnodes( innodesfile , odnodes );

        string outfile1 = outprefix+".ligand_target_ratio.txt";
        string outfile2 = outprefix+".sorted_ligand.txt";
        string outfile3 = outprefix+".target_ratio.txt";
        ana_ligand_target( bk, gs, odnodes, consmarker_gs, modulegene_gs, l_t_map, ligand_sender_receiver, outfile1, outfile2, outfile3, outprefix );


    } else{
        exit_with_help("Unkown function!");
    }

    return 1;
}

