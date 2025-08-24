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

} */


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

vector<string > parse_string_by_tab( string & instr )
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == '\t' )
		{
			
			strve.push_back(s);
			s = "";
			
		} else
		{
			s += instr[i];
		}
	}
	
	strve.push_back(s);
	
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

void transform_order( map<string, int > &tmap, multimap<int, string, greater<int> > &order_map )
{
    order_map.erase( order_map.begin(), order_map.end() );
    for ( map<string, int >::iterator ite = tmap.begin(); ite != tmap.end(); ++ite )
    {
        order_map.insert( make_pair( ite->second, ite->first ) );
    }
}

int sum_up( map<string, map<string, set<int > > > &tmp )
{
    int n = 0;
    for ( map<string, map<string, set<int > > >::iterator ite = tmp.begin(); ite != tmp.end(); ++ite )
    {
        for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            n += (int)si->second.size();
        }
    }
    return n;
}

class CPI
{
public:
    string id;
    string interacting_pair;
    string p_a;
    string p_b;
    string gene_a;
    string gene_b;
    bool secreted;
    bool receptor_a;
    bool receptor_b;
    bool is_integrin;
    string directionay;
    string classification;

};

class data_bank
{
public:
    vector<CPI> cpi_ve;
    map<string, int > ID_indexed_cpi;
    map<int, map<string, double > > ID_clusters_means;
    map<int, map<string, double > > ID_clusters_score;

    map<string, set<string> > complex_gene_map;
    map<string, string > gene_complex_map;

    map<string, set<int> > cl_cpiID;
    data_bank()
    {

    }

    void readin_sig_means_all(string infile );
    void readin_inter_score(string infile );
    void readin_deconv(string infile);
    void classify_cpi();

    map<string, map<string, set<string> > > site_cell_wingene;

    

};

void data_bank::readin_sig_means_all(string infile)
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    vector<string > tps = parse_string_by_tab(line);


    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string_by_tab(line);
        if ( ps.size() != tps.size() )
        {
            cout<<"error size "<<ps.size()<<" "<<tps.size()<<"\t"<<ps[0]<<endl; exit(1);
        }

        double rank = atof(ps[13].c_str());
        if ( rank > 1 )
        {
            cout<<"rank > 1"<<endl;
            break;
        }

        CPI cpi;
        cpi.id = ps[0];
        cpi.interacting_pair = ps[1];
        cpi.p_a = ps[2];
        cpi.p_b = ps[3];
        cpi.gene_a = ps[4];
        cpi.gene_b = ps[5];
        
        if ( ps[6] == "False")
            cpi.secreted = false;
        else if ( ps[6] == "True" )
            cpi.secreted = true;
        else
        {
            cout<<"error unexpected secreted "<<ps[6]<<" "<<cpi.id<<endl; exit(1);
        }
        if ( ps[7] == "False")
            cpi.receptor_a = false;
        else if ( ps[7] == "True")
            cpi.receptor_a = true;
        else
        {
            cout<<"error unexpected receptor_a "<<ps[7]<<" "<<cpi.id<<endl; exit(1);
        }
        if ( ps[8] == "False")
            cpi.receptor_b = false;
        else if ( ps[8] == "True")
            cpi.receptor_b = true;
        else
        {
            cout<<"error unexpected receptor_b "<<ps[8]<<" "<<cpi.id<<endl; exit(1);
        }
        if ( ps[10] == "False")
            cpi.is_integrin = false;
        else if ( ps[10] == "True")
            cpi.is_integrin = true;
        else
        {
            cout<<"error unexpected is_integrin "<<ps[10]<<" "<<cpi.id<<endl; exit(1);
        }
        cpi.directionay = ps[11];
        cpi.classification = ps[12];
        cpi_ve.push_back(cpi);
        int id = (int)cpi_ve.size()-1;

        ID_indexed_cpi.insert(make_pair(ps[0], id));

        for ( int i = 14; i < ps.size(); ++i )
        {
            string v = ps[i];
            if ( v != "" )
            {
                double means = atof(v.c_str() );
                string name = tps[i];
                if ( name.find("NP_HFC") == name.npos )
                    ID_clusters_means[id][name] = means;
            }
        }
        
    }
    inf.close();
}

void data_bank::readin_inter_score(string infile )
{

    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    cout<<"read file "<<infile<<endl;
    string line;
    getline(inf, line);
    vector<string > tps = parse_string_by_tab(line);

    while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string_by_tab(line);
        if ( ps.size() != tps.size() )
        {
            cout<<"error size "<<ps.size()<<" "<<tps.size()<<"\t"<<ps[0]<<endl; exit(1);
        }

        string id = ps[0];
        if ( ID_indexed_cpi.find(id) == ID_indexed_cpi.end( ) )
        {
            continue;
        }

        for ( int i = 13; i < ps.size(); ++i )
        {
            string v = ps[i];
            int idx = ID_indexed_cpi[id];
            double score = atof( ps[i].c_str() );
            if ( score > 0 )
            {
                if ( tps[i].find("NP_HFC") == tps[i].npos )
                    ID_clusters_score[idx][tps[i]] = score;
            }
        }


    }
    inf.close();
}

void data_bank::readin_deconv(string infile )
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
        gene_complex_map[gene] = comp;
        complex_gene_map[comp].insert(gene);
    }
    inf.close();

}

void data_bank::classify_cpi()
{
    for ( size_t i = 0; i < cpi_ve.size(); ++i )
    {
        string cl = cpi_ve[i].classification;
        if (cl =="")
            continue; 
        int idx = (int)i;
        cl_cpiID[cl].insert(idx);
    }
}

void readincolor(string infile, map<string, string > &color_map )
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
        vector<string > ps = parse_string(line);
        string cell = ps[0];
        string color = ps[1];
        color_map[cell] = color;
    }
    inf.close();
        
}

void sparse_name(string inname, string &site, string &cell_a, string &cell_b )
{

    vector<string > ps0 = parse_string(inname, '|' );
    vector<string > ps1_1 = parse_string( ps0[0], '_' );
    vector<string > ps1_2 = parse_string( ps0[1], '_' );
    site = ps1_1[0];
    cell_a = ps1_1[1];
    if ( (int)ps1_1.size() > 2 )
    {
        for ( size_t i = 2; i < ps1_1.size(); ++i )
            cell_a+=('_'+ps1_1[i]);
    }

    cell_b = ps1_2[1];
    if ( (int)ps1_2.size() > 2 )
    {
        for ( size_t i = 2; i < ps1_2.size(); ++i )
            cell_b+=('_'+ps1_2[i]);
    }

}

/*
void ana1( data_bank &db )
{
    ofstream outf("CCC_site.txt");
    
    map<string, map<string, map<string, vector<string > > > > site_a_b_cl;
    for ( map<int, map<string, double > >::iterator ite = db.ID_clusters_means.begin(); 
        ite != db.ID_clusters_means.end(); ++ite )
    {
        int idx = ite->first;
        string cl = db.cpi_ve[idx].classification;
        for ( map<string, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string name = si->first;
            string site = "";
            string cell_a = "";
            string cell_b = "";
            sparse_name(name, site, cell_a, cell_b );
            site_a_b_cl[site][cell_a][cell_b].push_back(cl);
        }
    }

    for ( map<string, map<string, map<string, vector<string > > > >::iterator ite = site_a_b_cl.begin(); 
        ite != site_a_b_cl.end(); ++ite )
    {
        string site = ite->first;
        outf<<">"<<site<<endl;
        for ( map<string, map<string, vector<string > > >::iterator si = ite->second.begin(); 
            si != ite->second.end(); ++si )
        {
            string cell_a = si->first;
            for ( map<string, vector<string > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
            {
                string cell_b = ti->first;
                outf<<cell_a<<"\t"<<cell_b<<endl;
                map<string, int > cl_c;
                
                for ( vector<string >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
                {
                    if ( cl_c.find(*fi) == cl_c.end() )
                        cl_c[*fi] = 1;
                    else
                        cl_c[*fi] += 1;
                }
                multimap<int, string, greater<int> > c_cl;

                for ( map<string, int >::iterator fi = cl_c.begin(); fi != cl_c.end(); ++fi )
                {
                    c_cl.insert(make_pair(fi->second, fi->first ) );
                }

                for ( multimap<int, string, greater<int> >::iterator fi = c_cl.begin(); fi != c_cl.end(); ++fi )
                {
                    outf<<fi->second<<":"<<fi->first<<endl;
                }
                
            }
        }
    }
    outf.close();
}

void ana2(data_bank &db)
{
    ofstream outf( "CCC_ana2.result.txt");
    for ( map<string, set<int> >::iterator ite = db.cl_cpiID.begin(); ite != db.cl_cpiID.end(); ++ite )
    {
        string cl = ite->first;
        outf<<">>"<<cl<<endl;
     //   cout<<cl<<endl;
        map<string, int > site_count;
        map<string, map<string, int> > site_cell_p_count;
        for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            int index = *si;
            if ( db.ID_clusters_means.find(index) == db.ID_clusters_means.end() )
            {
                cout<<"error cannot find index "<<index<<endl; exit(1);
            }
            for ( map<string, double >::iterator pi = db.ID_clusters_means[index].begin(); 
                pi != db.ID_clusters_means[index].end(); ++pi )
            {
                string name = pi->first;
                string site = "";
                string cell_a = "";
                string cell_b = "";
                sparse_name(name, site, cell_a, cell_b );
                if ( site_count.find(site) == site_count.end() )
                {   
                    site_count[site] = 1;
                    site_cell_p_count[site][name] = 1;
                } else
                {
                    site_count[site] += 1;
                    if ( site_cell_p_count[site].find( name ) == site_cell_p_count[site].end() )
                    {
                        site_cell_p_count[site][name] = 1;
                    } else{
                        site_cell_p_count[site][name] += 1;
                    }
                }


            }
        } 
        
        for ( map<string, int >::iterator ssi = site_count.begin(); ssi != site_count.end(); ++ssi )
        {
            if ( ssi == site_count.begin() )
            {
                outf<<">";
            } else
            {
                outf<<"\t";
            }
            outf<<ssi->first<<":"<<ssi->second;
        }
        outf<<endl;

        for ( map<string, map<string, int> >::iterator ssi = site_cell_p_count.begin(); 
            ssi != site_cell_p_count.end(); ++ssi )
        {
            outf<<"<"<<ssi->first;
            multimap<int, string, greater<int> > c_cellp;
            for ( map<string, int>::iterator ti = ssi->second.begin(); ti != ssi->second.end(); ++ti )
            {
                c_cellp.insert( make_pair(ti->second, ti->first ) );
            }
            for ( multimap<int, string, greater<int> >::iterator ti = c_cellp.begin(); ti != c_cellp.end(); ++ti )
            {
                outf<<"\t"<<ti->second<<":"<<ti->first;
            }
            outf<<endl;
        }
    }
    outf.close();
}

void ana3(data_bank &db, map<string, string > &color_map )
{
    map<string, map<string, set<int > > > cell_pair_inter;
    map<string, map<string, map<string, set<int > > > > site_cell_pair_inter;
    for ( map<int, map<string, double > >::iterator ite = db.ID_clusters_means.begin(); ite != db.ID_clusters_means.end(); ++ite )
    {
        int index = ite->first;
        for ( map<string, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string name = si->first;
            string site = "";
            string c_a = "";
            string c_b = "";
            sparse_name( name, site, c_a, c_b );
            cell_pair_inter[c_a][c_b].insert(index);
            site_cell_pair_inter[site][c_a][c_b].insert(index);
        }
    }

    string outfile1 = "total_cell_pair_ccc.df.txt";
    string outfile2 = "total_ligand.cl.txt";
    string outfile3 = "total_receptor.cl.txt";
    ofstream outf1(outfile1.data() );
    outf1<<"ligend\treceptor\tinteractions"<<endl;

    map<string, set<int> > l_c;
    map<string, set<int> > r_c;
    map<string, int > l_count;
    map<string, int > r_count;
    for ( map<string, map<string, set<int > > >::iterator ite = cell_pair_inter.begin(); ite != cell_pair_inter.end(); ++ite )
    {
        string c_a = ite->first;
        for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string c_b = si->first;
            int n = (int)si->second.size();
            outf1<<c_a<<"_l\t"<<c_b<<"_r\t"<<n<<endl;
            l_c[c_a].insert(si->second.begin(), si->second.end());
            r_c[c_b].insert(si->second.begin(), si->second.end());

            if ( l_count.find(c_a) == l_count.end() )
            {
                l_count[c_a] = n;
            } else
            {
                l_count[c_a] += n;
            }
            if ( r_count.find(c_b) == r_count.end() )
            {
                r_count[c_b] = n;
            } else
            {
                r_count[c_b] += n;
            }
        }
    }
    outf1.close();

    multimap<int, string, greater<int> > order_l;
    multimap<int, string, greater<int > > order_r;
    for ( map<string, int >::iterator ite = l_count.begin(); ite != l_count.end(); ++ite )
    {
        order_l.insert( make_pair(ite->second, ite->first ));
    }
    for ( map<string, int >::iterator ite = r_count.begin(); ite != r_count.end(); ++ite )
    {
        order_r.insert( make_pair(ite->second, ite->first ));
    }
    ofstream outf_orderl( "total_ligend_order.txt");
    ofstream outf_orderr( "total_receptor_order.txt");
    for ( multimap<int, string, greater<int> >::iterator ite = order_l.begin(); ite != order_l.end(); ++ite )
    {
        outf_orderl<<ite->second<<"_l"<<endl;
    }
    for ( multimap<int, string, greater<int> >::iterator ite = order_r.begin(); ite != order_r.end(); ++ite )
    {
        outf_orderr<<ite->second<<"_r"<<endl;
    }
    outf_orderl.close();
    outf_orderr.close();

    ofstream outf2(outfile2.data() );
    for ( map<string, set<int> >::iterator ite = l_c.begin(); ite != l_c.end(); ++ite )
    {
        string c = ite->first;
        map<string, int > cl_c;
        for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            int index = *si;
            string cl = db.cpi_ve[index].classification;
            if ( cl_c.find(cl) == cl_c.end() )
            {
                cl_c[cl] = 1;
            } else
            {
                cl_c[cl] += 1;
            }
        }
        multimap<int, string, greater<int> > order_cl;
        for ( map<string, int >::iterator si = cl_c.begin(); si != cl_c.end(); ++si )
        {
            order_cl.insert( make_pair(si->second, si->first ) );
        }
        outf2<<">>"<<c<<endl;
        for ( multimap<int, string, greater<int> >::iterator si = order_cl.begin(); si != order_cl.end(); ++si )
        {
            outf2<<si->second<<"\t"<<si->first<<endl;
        }

    }
    outf2.close();

    ofstream outf3(outfile3.data() );
    for ( map<string, set<int> >::iterator ite = r_c.begin(); ite != r_c.end(); ++ite )
    {
        string c = ite->first;
        map<string, int > cl_c;
        for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            int index = *si;
            string cl = db.cpi_ve[index].classification;
            if ( cl_c.find(cl) == cl_c.end() )
            {
                cl_c[cl] = 1;
            } else
            {
                cl_c[cl] += 1;
            }
        }
        multimap<int, string, greater<int> > order_cl;
        for ( map<string, int >::iterator si = cl_c.begin(); si != cl_c.end(); ++si )
        {
            order_cl.insert( make_pair(si->second, si->first ) );
        }
        outf3<<">>"<<c<<endl;
        for ( multimap<int, string, greater<int> >::iterator si = order_cl.begin(); si != order_cl.end(); ++si )
        {
            outf3<<si->second<<"\t"<<si->first<<endl;
        }

    }
    outf3.close();

    map<string, map<string, int > > cell_pair_count;
    int total_count = 0;
    map<string, int > site_count;
    
    for ( map<string, map<string, map<string, set<int > > > >::iterator ite = site_cell_pair_inter.begin();
        ite != site_cell_pair_inter.end(); ++ite )
    {
        int c = sum_up(ite->second);
        site_count[ite->first] = c;
        total_count += c;
        int exclude_n = 0;
        for ( map<string, map<string, set<int > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string c_a = si->first;
            for (map<string, set<int > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
            {
                string c_b = ti->first;
                int n = (int)ti->second.size();
                if ( c_a == "HFC" | c_b == "HFC" | c_a == "Sebocytes" | c_b == "Sebocytes" )
                {
                    exclude_n += n;
                } 

                if ( cell_pair_count[c_a].find( c_b ) == cell_pair_count[c_a].end() )
                {
                    cell_pair_count[c_a][c_b] = n;
                } else
                {
                    cell_pair_count[c_a][c_b] += n;
                }
                
            }
        }

        total_count -= exclude_n;
        site_count[ite->first] -= exclude_n;
    }

    for ( map<string, map<string, map<string, set<int > > > >::iterator ite = site_cell_pair_inter.begin();
        ite != site_cell_pair_inter.end(); ++ite )
    {
        string site = ite->first;
        string outfile_s = site+"_cell_pair_ccc.df.txt";
        ofstream outf_s( outfile_s.data() );

        map<string, int > l_count_s;
        map<string, int > r_count_s;  

        outf_s<<"ligend\treceptor\tinteractions\tHyperGp\tColor"<<endl;
        for ( map<string, map<string, set<int > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string c_a = si->first;
            for ( map<string, set<int > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
            {
                string c_b = ti->first;
                int c= (int)ti->second.size();
                int r = cell_pair_count[c_a][c_b];
                int n = site_count[site];
                
                double hyperg = hypergeometrictest(r, n, total_count, c);
                if ( c_a == "HFC" | c_b == "HFC" | c_a == "Sebocytes" | c_b == "Sebocytes" )
                {
                    hyperg = 1;
                }

                string color = color_map[c_a];
                if ( hyperg < 0.05 )
                    color += "FF";
                else
                    color += "30";

                outf_s<<c_a<<"_l\t"<<c_b<<"_r\t"<<c<<"\t"<<hyperg<<"\t"<<color<<endl;
                

                if ( l_count_s.find(c_a) == l_count_s.end())
                {
                    l_count_s[c_a] = c;
                } else
                {
                    l_count_s[c_a] += c;
                }
                if ( r_count_s.find(c_b) == r_count_s.end() )
                {
                    r_count_s[c_b] = c;
                } else
                {
                    r_count_s[c_b] += c;
                }
            }
        }
        outf_s.close();
        
        multimap<int, string, greater<int> > order_l_s;
        multimap<int, string, greater<int> > order_r_s;
        transform_order( l_count_s, order_l_s );
        transform_order( r_count_s, order_r_s );
        string outfile_s2 = site+"_ligend_order.txt";
        ofstream outf_s2(outfile_s2.data( ));
        for ( multimap<int, string, greater<int> >::iterator si = order_l_s.begin(); si != order_l_s.end(); ++si )
        {
            outf_s2<<si->second<<"_l"<<endl;
        }
        outf_s2.close();
        string outfile_s3 = site+"_receptor_order.txt";
        ofstream outf_s3(outfile_s3.data( ));
        for ( multimap<int, string, greater<int> >::iterator si = order_r_s.begin(); si != order_r_s.end(); ++si )
        {
            outf_s3<<si->second<<"_r"<<endl;
        }
        outf_s3.close();
    }

}
*/
void readin_DEgene(string infile, string celltype, map<string, map<string, set<string> > > &site_cell_wingene )
{
    ifstream inf(infile.data() );

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
        vector<string > ps = parse_string(line);
        string gene = ps[0];
        string site = ps[1];
        site_cell_wingene[site][celltype].insert(gene);
    }
    inf.close();
}

void readin_DEgene_b(string infile, map<string, map<string, set<string> > > &site_cell_wingene )
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
        vector<string > ps = parse_string(line, ' ');
        string cell = ps[0];
        string ifile = ps[1]; 
        readin_DEgene( ifile, cell, site_cell_wingene );
    }
    inf.close();
}

void readinmodule( string infile, string name, map<string, map<string, string > > &genemodules ) // celltype, gene, module
{

	ifstream inf(infile.data() );

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
		vector<string > ps = parse_string( line, ',' );
		string gene = ps[0];
		string m = ps[5];
     //   if ( m != "0" )
     //   {
            m = "Module"+m;
            genemodules[name][gene] = m;
     //   }
		
	}
	inf.close();

}

void readinmodule_b( string infile, map<string, map<string, string > > &genemodules )
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
		vector<string > ps = parse_string( line, '_' );
        string name = ps[1];
        readinmodule( line, name, genemodules );
    }
    inf.close();

	
}

void ana4( data_bank &db, map<string, map<string, string > > &genemodules, string outfile )  // intersecting ccc with module genes.
{

    map<string, map<string, string > > gene_celltype_modules;
    for ( map<string, map<string, string > >::iterator ite = genemodules.begin(); 
        ite != genemodules.end(); ++ite )
    {
        string celltype = ite->first;
        for ( map<string, string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string gene = si->first;
            string m = si->second;
            gene_celltype_modules[gene][celltype] = m;
        }
    }

    
    ofstream outf( outfile.data() );
    outf<<"Cid\tGene_a\tCelltype_a\tModule_a\tGene_b\tCelltype_b\tModule_b"<<endl;

    for ( size_t i = 0; i < db.cpi_ve.size(); ++i )
    {

        int id = (int)i;
        map<string, double > name_mean = db.ID_clusters_means[i];
        map<string, set<string > > celltypes_a_parts;
        map<string, set<string > > celltypes_b_parts;
        for ( map<string, double >::iterator ite = name_mean.begin(); ite != name_mean.end(); ++ite )
        {
            string name = ite->first;
            string site = "";
            string cell_a = "";
            string cell_b = "";
            sparse_name(name, site, cell_a, cell_b );
            celltypes_a_parts[cell_a].insert( cell_b );
            celltypes_b_parts[cell_b].insert( cell_a );
        }

        map<pair<string, string >, set<pair<string, string > >  > paired_genes;  // <celltype, gene> <celltype, gene>

        set<string > geneset_a;
        set<string > geneset_b;
        if ( db.cpi_ve[i].gene_a != "" )
        {
            geneset_a.insert(db.cpi_ve[i].gene_a);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_a, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_a<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_a.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_a = db.complex_gene_map[comp];
                }
            }

        }

        if ( db.cpi_ve[i].gene_b != "" )
        {
            geneset_b.insert(db.cpi_ve[i].gene_b);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_b, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_b<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_b.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_b = db.complex_gene_map[comp];
                }
            }
        }

        for ( set<string >::iterator ite_a = geneset_a.begin(); ite_a != geneset_a.end(); ++ite_a )
        {
            string ga = *ite_a;
            if ( gene_celltype_modules.find( ga ) == gene_celltype_modules.end() )
            {
                continue;
            } 
            for ( set<string>::iterator ite_b = geneset_b.begin(); ite_b != geneset_b.end(); ++ite_b )
            {
                string gb = *ite_b;
                if ( gene_celltype_modules.find( gb ) == gene_celltype_modules.end() )
                    continue;

                for ( map<string, set<string > >::iterator ci = celltypes_a_parts.begin(); 
                    ci != celltypes_a_parts.end(); ++ci )
                {
                    string ca = ci->first;
                    if ( gene_celltype_modules[ga].find( ca ) == gene_celltype_modules[ga].end() )
                        continue;
                    for ( set<string >::iterator sci = ci->second.begin(); sci != ci->second.end(); ++sci )
                    {
                        string cb = *sci;
                        if ( gene_celltype_modules[gb].find( cb ) == gene_celltype_modules[gb].end() )
                        {
                            continue;
                        } 
                        paired_genes[make_pair(ca, ga)].insert( make_pair(cb, gb) );
                    }
                }
            }

            if ( !paired_genes.empty() )
            {
                for ( map<pair<string, string >, set<pair<string, string > >  >::iterator pi = paired_genes.begin();
                    pi != paired_genes.end(); ++pi )
                {
                    string ca = pi->first.first;
                    string ga = pi->first.second;
                    string ma = gene_celltype_modules[ga][ca];
                    for ( set<pair<string, string > >::iterator spi = pi->second.begin(); spi != pi->second.end(); ++spi )
                    {
                        string cb = spi->first;
                        string gb = spi->second;
                        string mb = gene_celltype_modules[gb][cb];
                        outf<<db.cpi_ve[i].id<<"\t"<<ca<<"\t"<<ga<<"\t"<<ma<<"\t"<<cb<<"\t"<<gb<<"\t"<<mb<<endl;
                    }
                }
                
            }
        }
    }
    outf.close();
}

void readin_DEgene_2( string infile, map<string, set<string > > &gene_celltypes)
{
    ifstream inf(infile.data() );

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
        vector<string > ps = parse_string(line);
        string gene = ps[1];
        string ct = ps[3];
        vector<string > ct_ve = parse_string( ct, ',');
        for ( size_t i = 0; i < ct_ve.size(); ++i )
        {
            string celltype = ct_ve[i];
            gene_celltypes[gene].insert(celltype);
        }
    }
    inf.close();
}

void ana5( data_bank &db, map<string, set<string > > &gene_celltypes, string outfile )  // intersecting ccc with module genes.
{

    
    ofstream outf( outfile.data() );
    outf<<"Cid\tGene_a\tCelltype_a\tGene_b\tCelltype_b"<<endl;

    for ( size_t i = 0; i < db.cpi_ve.size(); ++i )
    {

        int id = (int)i;
        map<string, double > name_mean = db.ID_clusters_means[i];
        map<string, set<string > > celltypes_a_parts;
        map<string, set<string > > celltypes_b_parts;
        for ( map<string, double >::iterator ite = name_mean.begin(); ite != name_mean.end(); ++ite )
        {
            string name = ite->first;
            string site = "";
            string cell_a = "";
            string cell_b = "";
            sparse_name(name, site, cell_a, cell_b );
            celltypes_a_parts[cell_a].insert( cell_b );
            celltypes_b_parts[cell_b].insert( cell_a );
        }

        map<pair<string, string >, set<pair<string, string > >  > paired_genes;  // <celltype, gene> <celltype, gene>

        set<string > geneset_a;
        set<string > geneset_b;
        if ( db.cpi_ve[i].gene_a != "" )
        {
            geneset_a.insert(db.cpi_ve[i].gene_a);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_a, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_a<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_a.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_a = db.complex_gene_map[comp];
                }
            }

        }

        if ( db.cpi_ve[i].gene_b != "" )
        {
            geneset_b.insert(db.cpi_ve[i].gene_b);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_b, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_b<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_b.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_b = db.complex_gene_map[comp];
                }
            }
        }

        for ( set<string >::iterator ite_a = geneset_a.begin(); ite_a != geneset_a.end(); ++ite_a )
        {
            string ga = *ite_a;
            if ( gene_celltypes.find( ga ) == gene_celltypes.end() )
            {
                continue;
            } 
            for ( set<string>::iterator ite_b = geneset_b.begin(); ite_b != geneset_b.end(); ++ite_b )
            {
                string gb = *ite_b;
                if ( gene_celltypes.find( gb ) == gene_celltypes.end() )
                    continue;

                for ( map<string, set<string > >::iterator ci = celltypes_a_parts.begin(); 
                    ci != celltypes_a_parts.end(); ++ci )
                {
                    string ca = ci->first;
                    if ( gene_celltypes[ga].find( ca ) == gene_celltypes[ga].end() )
                        continue;
                    for ( set<string >::iterator sci = ci->second.begin(); sci != ci->second.end(); ++sci )
                    {
                        string cb = *sci;
                        if ( gene_celltypes[gb].find( cb ) == gene_celltypes[gb].end() )
                        {
                            continue;
                        } 
                        paired_genes[make_pair(ca, ga)].insert( make_pair(cb, gb) );
                    }
                }
            }

            if ( !paired_genes.empty() )
            {
                for ( map<pair<string, string >, set<pair<string, string > >  >::iterator pi = paired_genes.begin();
                    pi != paired_genes.end(); ++pi )
                {
                    string ca = pi->first.first;
                    string ga = pi->first.second;
                    for ( set<pair<string, string > >::iterator spi = pi->second.begin(); spi != pi->second.end(); ++spi )
                    {
                        string cb = spi->first;
                        string gb = spi->second;
                        outf<<db.cpi_ve[i].id<<"\t"<<ca<<"\t"<<ga<<"\t"<<cb<<"\t"<<gb<<"\t"<<endl;
                    }
                }
                
            }
        }
    }
    outf.close();
}

string rename(string inname )
{
    if ( inname == "KC_Channel" )
        return "Channel";
    else if ( inname == "Pc-vSMC" )
        return "Pc";
    else if ( inname == "Mac-DC" )
        return "Mac";
    else 
        return inname;
}

void ana6( data_bank &db, string outfile )  // intersecting ccc with module genes.
{

    
    ofstream outf( outfile.data() );
    outf<<"Gene_a\tCelltype_a\tGene_b\tCelltype_b"<<endl;

    map<pair<string, string >, set<pair<string, string > >  > paired_genes;  // <celltype, gene> <celltype, gene>

    for ( size_t i = 0; i < db.cpi_ve.size(); ++i )
    {

        int id = (int)i;
        map<string, double > name_mean = db.ID_clusters_means[i];
        map<string, set<string > > celltypes_a_parts;
        map<string, set<string > > celltypes_b_parts;
        for ( map<string, double >::iterator ite = name_mean.begin(); ite != name_mean.end(); ++ite )
        {
            string name = ite->first;
            string site = "";
            string cell_a = "";
            string cell_b = "";
            sparse_name(name, site, cell_a, cell_b );
            celltypes_a_parts[cell_a].insert( cell_b );
            celltypes_b_parts[cell_b].insert( cell_a );
        }

        

        set<string > geneset_a;
        set<string > geneset_b;
        if ( db.cpi_ve[i].gene_a != "" )
        {
            geneset_a.insert(db.cpi_ve[i].gene_a);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_a, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_a<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_a.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_a = db.complex_gene_map[comp];
                }
            }

        }

        if ( db.cpi_ve[i].gene_b != "" )
        {
            geneset_b.insert(db.cpi_ve[i].gene_b);
        } else
        {
            vector<string > ps = parse_string(db.cpi_ve[i].p_b, ':' );
            if ( ps[0] != "complex" )
            {
                cout<<"unexpected gene and p "<<db.cpi_ve[i].p_b<<endl;
                exit(1);
            } else
            {
                string comp = db.cpi_ve[i].p_b.substr(8);
                if ( db.complex_gene_map.find(comp) == db.complex_gene_map.end() )
                {
                    cout<<"error cannot find comp in complexgenemap "<<comp<<endl; exit(1);
                } else
                {
                    geneset_b = db.complex_gene_map[comp];
                }
            }
        }

        for ( set<string >::iterator ite_a = geneset_a.begin(); ite_a != geneset_a.end(); ++ite_a )
        {
            string ga = *ite_a;
            
            for ( set<string>::iterator ite_b = geneset_b.begin(); ite_b != geneset_b.end(); ++ite_b )
            {
                string gb = *ite_b;
                

                for ( map<string, set<string > >::iterator ci = celltypes_a_parts.begin(); 
                    ci != celltypes_a_parts.end(); ++ci )
                {
                    string ca = ci->first;
                    
                    for ( set<string >::iterator sci = ci->second.begin(); sci != ci->second.end(); ++sci )
                    {
                        string cb = *sci;
                        
                        paired_genes[make_pair(ca, ga)].insert( make_pair(cb, gb) );

                        if ( ca == "Fb" && ga == "TENM3" && cb == "Fb" && gb == "ADGRL3")
                        {
                            cout<<"echo "<<db.cpi_ve[i].id<<endl; 
                        }
                    }
                }
            }

            
        }
    }

    if ( !paired_genes.empty() )
    {
        for ( map<pair<string, string >, set<pair<string, string > >  >::iterator pi = paired_genes.begin();
            pi != paired_genes.end(); ++pi )
        {
            string ca = pi->first.first;
            string ga = pi->first.second;
            for ( set<pair<string, string > >::iterator spi = pi->second.begin(); spi != pi->second.end(); ++spi )
            {
                string cb = spi->first;
                string gb = spi->second;
                outf<<rename(ca)<<"\t"<<ga<<"\t"<<rename(cb)<<"\t"<<gb<<"\t"<<endl;
            }
        }
        
    }
    outf.close();
}

int main(int argc, char* argv[])
{
    if ( argc == 1 )
    {
        cout<<"analyze cellphone db results "<<endl;
//        cout<<"Usage: prog sig_means scores deconvfile integratedDEgene outfile"<<endl;
        cout<<"Usage: prog sig_means scores deconvfile outfile"<<endl;

        exit(1);
    }

    string infile1 = argv[1];
    string infile2 = argv[2];
 //   string colorfile = argv[3];
    string infile3 = argv[3];
 //   string genemodulefile = argv[4];
 //   string degenefile = argv[4];
    string outfile = argv[4];
    data_bank db;

    db.readin_sig_means_all( infile1 );
    db.readin_inter_score( infile2 );
    db.readin_deconv(infile3);
    db.classify_cpi();

 //   ana1(db);
 //   
 //   ana2(db);

 //   map<string, string > color_map;
 //   readincolor( colorfile, color_map );
 //   ana3(db, color_map);

 //   map<string, map<string, string > > genemodules;
 //   readinmodule_b( genemodulefile, genemodules );
 //   ana4( db, genemodules, outfile );

 /*   map<string, set<string > > gene_celltypes;
    readin_DEgene_2( degenefile, gene_celltypes );

    ana5( db, gene_celltypes, outfile );  */

    ana6( db, outfile );

    return 1;
}
