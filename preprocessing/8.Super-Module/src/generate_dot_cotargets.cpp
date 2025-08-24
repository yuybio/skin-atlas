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

using namespace std;

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

//    vector<CellChat> CC_ve;

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
//    void readincellchat( string infile );
//    void readincellchat3( string infile );
//    void readindeconvolute( string infile );
//    void readinoverlappedgenes( string infile );

 //   void deconvolutegetset();
//    void AssignCellchat2link( Bank &bk );
//    void AssignCellchat2link2( Bank &bk );

//    void readinsigovlgenes( string infile, double pthr, int type );
//    void readinsharedgokegg( string infile, set<string > &proi_go );
//    void readincellchat2( string infile );

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

void readinsharedtargetfile( string infile, map<string, DOT_g > &ligand_dots )
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
        string nodea = ps[0];
        string nodeb = ps[1];
        string cotargets = ps[9];
        if ( cotargets == "NA" )
            continue;
        vector<string> ps_1 = parse_string( cotargets, '|');
        for ( size_t i = 0; i < ps_1.size(); ++i )
        {
            string cotargets_sp1 = ps_1[i];
            vector<string > ps_2 = parse_string( cotargets_sp1, ':');
            string ligand = ps_2[0];
            string cotargets_sp2 = ps_2[1];
            vector<string > targets_ve = parse_string( cotargets_sp2, ',');
            if ( ligand_dots.find( ligand ) != ligand_dots.end() )
            {
                for ( size_t j = 0; j < targets_ve.size(); ++j )
                {
                    string t = targets_ve[j];
                    ligand_dots[ligand].targets.insert(t);
                    ligand_dots[ligand].ligand_target[ligand].insert(t);
                    ligand_dots[ligand].target_node[t].insert(nodea);
                    ligand_dots[ligand].target_node[t].insert(nodeb);
                }
            } else{
                DOT_g dg;
                dg.ligands.insert( ligand );
                for ( size_t j = 0; j < targets_ve.size(); ++j )
                {
                    string t = targets_ve[j];
                    dg.targets.insert(t);
                    dg.ligand_target[ligand].insert(t);
                    dg.target_node[t].insert(nodea);
                    dg.target_node[t].insert(nodeb);
                }
                ligand_dots.insert( make_pair(ligand, dg) );
            }
        }
    }
    inf.close();
}

void readinligandnode( string infile, map<string, set<string > > &ligand_nodes) 
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
        char firstc = line[0];
        if ( firstc != '>' )
            continue;
        vector<string > ps = parse_string( line, '\t' );
        string ligand = ps[0].substr(1);
        string node_s = ps[1];
        vector<string > node_ve = parse_string( node_s, ':');
        for ( size_t i = 0; i < node_ve.size(); ++i )
        {
            ligand_nodes[ligand].insert( node_ve[i]);
        }
    }
    inf.close();
}

void addinligandnode( map<string, DOT_g > &ligand_dots, map<string, set<string > > &ligand_nodes )
{
    for ( map<string, DOT_g >::iterator ite = ligand_dots.begin(); ite != ligand_dots.end(); ++ite )
    {
        string ligand = ite->first;
        if ( ligand_nodes.find(ligand ) == ligand_nodes.end() )
        {
            cout<<"error cannot find nodes for ligand "<<ligand<<endl; exit(1);
        }
        ite->second.ligand_node[ligand].insert( ligand_nodes[ligand].begin(), ligand_nodes[ligand].end() );
    }
}


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
   
    if ( n_o >= thr_n )
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

void integrate(  string prefix, Bank &bk, GeneSet &gs, map<string, DOT_g > &ligand_dots )
{
    vector<DOT_g > dot_ve;
    for ( map<string, DOT_g >::iterator ite = ligand_dots.begin(); ite != ligand_dots.end(); ++ite )
    {
        string ligand = ite->first;
        dot_ve.push_back( ite->second );
    }

    int ovl_n_thr = 1;
    vector<DOT_g> merged_DOT_g;
    merge_DOT_g( dot_ve, merged_DOT_g, ovl_n_thr );

    cout<<"generate_dotfile"<<endl;
    generate_dotfile( merged_DOT_g, prefix, bk, gs );
}



int main(int argc, char* argv[] )
{
    if ( argc == 1 )
    {
        cout<<"generate dot plot for co-targets"<<endl; 
        cout<<"Usage: prog infile1[conservedmarkers] infile2[allmodules] infile3[node] infile4[link] infile5[co-targetfile] infile6[sorted-ligandfile] outprefix"<<endl;
        exit(1);
    }

    string inconsgenefile = argv[1];
    string inmodulegenefile = argv[2];
    string innodesfile = argv[3];
    string inlinksfile = argv[4];

    Bank bk;
    bk.readinconsgeneset( inconsgenefile );

    bk.readinmodulegeneset( inmodulegenefile );

    GeneSet gs;
    gs.readingeneset( innodesfile );

    gs.readinlinks( inlinksfile );

    string incotargetfile = argv[5];
    map<string, DOT_g > ligand_dots;
    readinsharedtargetfile( incotargetfile, ligand_dots );

    string inligandnodefile = argv[6];
    map<string, set<string > > ligand_nodes;
    readinligandnode(inligandnodefile, ligand_nodes);

    addinligandnode( ligand_dots, ligand_nodes );

    string prefix = argv[7];
    integrate( prefix, bk, gs, ligand_dots );

    return 1;
}
