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

void readin( string infile, map<string, map<string, double > > &mod_samp_exp )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;

    size_t p = infile.find("_module");
    string name = infile.substr(0, (int)p);
    cout<<name<<endl;
    string line;
	getline(inf, line);
    vector<string > mods = parse_string( line, ',' );
    
    
    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line, ',' );
        string samp = ps[0];
        if ( ps.size() != mods.size()+1 )
        {
            cout<<"error ps "<<ps.size()<<" "<<mods.size()<<endl; exit(1);
        }
        for ( size_t i = 1; i < ps.size(); ++i )
        {
            double sc = atof(ps[i].c_str() );
            string m = name+"_"+mods[i-1];
            mod_samp_exp[m][samp] = sc;
           
        }
    }
    inf.close();
   
}

void readinbatch( string infile, map<string, map<string, double > > &mod_samp_exp )
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
        string filename = line;
        readin( filename, mod_samp_exp );
    }
    inf.close();

}

void matrix_output( string outfile, map<string, map<string, double > > &mod_samp_exp )
{
    ofstream outf( outfile.data() );
    outf<<"Mod";
    set<string > samples;
    for ( map<string, map<string, double > >::iterator ite = mod_samp_exp.begin(); 
        ite != mod_samp_exp.end(); ++ite )
    {
    //    cout<<ite->first<<endl;
        for ( map<string, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            outf<<"\t"<<si->first;
            samples.insert(si->first);
         //   cout<<si->first<<endl;
        }
        outf<<endl;
        break;

    }
    

    for ( map<string, map<string, double > >::iterator ite = mod_samp_exp.begin(); 
        ite != mod_samp_exp.end(); ++ite )
    {
        string mod = ite->first;
        outf<<mod;
        for ( set<string>::iterator si = samples.begin(); si != samples.end(); ++si )
        {
            double v = 0;
            if ( ite->second.find( *si ) != ite->second.end() )
            {
                v = ite->second[*si];
                
            }
            outf<<"\t"<<v;
        }
        outf<<endl;
    }
    outf.close();
}

int main( int argc, char* argv[] )
{
    if ( argc == 1 )
    {
        cout<<"combine all module to one matrix for correlation"<<endl;
        cout<<"Usage: prog infile outfile "<<endl;
        exit(1);
    }

    string infile = argv[1];
    string outfile = argv[2];

    map<string, map<string, double > > mod_samp_exp;
    readinbatch( infile, mod_samp_exp );

    matrix_output( outfile, mod_samp_exp );

    return 1;
}

