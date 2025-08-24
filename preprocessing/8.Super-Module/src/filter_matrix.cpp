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

void readinMatrix( string infile, map<string, map<string, double > > &cor_matrix )
{
	
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
	getline(inf, line);
	vector<string > name_a = parse_string( line, '\t' );


    
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, '\t' );
		string nb = ps[0];
		for ( size_t i = 1; i < ps.size(); ++i )
		{
			double cor = atof(ps[i].c_str() );
			string na = name_a[i-1];
			cor_matrix[na][nb] = cor;
			cor_matrix[nb][na] = cor;
			
			
		}

	}
	inf.close();
}

void readinterms( string infile, set<string > &terms )
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
        vector<string > ps = parse_string( line, '\t' );
		string t = ps[0];
		terms.insert(t);
	}
	inf.close();
}

void filterandoutput( map<string, map<string, double > > &cor_matrix,
	map<string, map<string, double > > &p_matrix,
	set<string > &terms, 
	double cor_thr,
	double p_thr,
	string outfile1,
	string outfile2 )
{
    ofstream outf1(outfile1.data());
    ofstream outf2(outfile2.data());
    outf1<<"source\ttarget\tcor\tp"<<endl;
    map<string, double> agg_cor;
	for ( set<string >::iterator ite = terms.begin(); ite != terms.end(); ++ite )
	{
		string ta = *ite;
		if ( cor_matrix.find(ta) == cor_matrix.end() )
		{
			cout<<"error cannot find "<<ta<<" in cor_matrix"<<endl; exit(1);
		}
		set<string >::iterator si = ite;
		++si;
		for ( ; si != terms.end(); ++si )
		{
			string tb = *si;
			if ( cor_matrix[ta].find(tb) == cor_matrix[ta].end() )
			{
				cout<<"error cannot find "<<tb<<" in cor_matrix "<<ta<<endl; exit(1);
			}
			double cor = cor_matrix[ta][tb];
			double p = p_matrix[ta][tb];
			
            if ( cor >= cor_thr && p <= p_thr )
            {
                outf1<<ta<<"\t"<<tb<<"\t"<<cor<<"\t"<<p<<endl;
                if ( agg_cor.find(ta) != agg_cor.end() )
                {
                    agg_cor[ta] += cor;
                } else
                {
                    agg_cor[ta] = cor;
                }
                if ( agg_cor.find(tb) != agg_cor.end() )
                {
                    agg_cor[tb] += cor;
                } else{
                    agg_cor[tb] = cor;
                }
            }
		}
		
	}
    outf1.close();
    outf2<<"name\taggsize"<<endl;
    for ( map<string, double >::iterator ite = agg_cor.begin(); ite != agg_cor.end(); ++ite )
    {
        outf2<<ite->first<<"\t"<<ite->second<<endl;
    }
    outf2.close();
}

int main(int argc, char* argv[] )
{
    if ( argc != 8 )
    {
        cout<<"transform to network plot with filtering sig edge"<<endl;
        cout<<"Usage: prog in_cor_matrix in_p_matrix inselectterm cor_thr[double] p_thr[double] outfile1[link] outfile2[vert]"<<endl;
        exit(1);
    }

    string incorfile = argv[1];
    string inpfile = argv[2];
    string intermfile = argv[3];
    double cor_thr = atof(argv[4]);
    double p_thr = atof(argv[5]);
    string outlinkfile = argv[6];
    string outvertfile = argv[7];

	cout<<cor_thr<<" "<<p_thr<<endl;
    map<string, map<string, double > > cor_matrix;
    readinMatrix( incorfile, cor_matrix );

    map<string, map<string, double > > p_matrix;
    readinMatrix( inpfile, p_matrix );

    set<string > terms;
    readinterms( intermfile, terms );

    filterandoutput( cor_matrix, p_matrix, terms, cor_thr, p_thr, outlinkfile, outvertfile );

}

