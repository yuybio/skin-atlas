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

void readinlinkcor( string infile, map<string, map<string, pair<double,double> > > &cor_matrix )
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
        vector<string > ps = parse_string( line, '\t' );
		string na = ps[0];
		string nb = ps[1];
		double cor = atof(ps[2].c_str() );
        double corp = atof(ps[3].c_str() );
		cor_matrix[na][nb] = make_pair(cor,corp);
		cor_matrix[nb][na] = make_pair(cor,corp);

	}
	inf.close();
}

void readinvert( string infile, map<string, string > &node_celltype )
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
        vector<string > ps = parse_string( line, '\t' );
        string node = ps[0];
        string celltype = ps[2];
        node_celltype[node]= celltype;

    }
    inf.close();
}

void seperate( map<string, map<string, pair<double,double> > > &cor_matrix, 
    map<string, string > &node_celltype,
    map<string, set<string > > &sig_pairs,
    map<string, set<string > > &nonsig_pairs,
    double pthr, double corthr )
{
    for ( map<string, map<string, pair<double,double> > >::iterator ite = cor_matrix.begin(); ite != cor_matrix.end(); ++ite )
    {
        string nodea = ite->first;
        string cta = node_celltype[nodea];
        for ( map<string, pair<double, double > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string nodeb = si->first;
            string ctb = node_celltype[nodeb];
            if ( cta == ctb )
                continue;
            double cor = si->second.first;
            double p = si->second.second;
            if ( cor >= corthr && p <= pthr )
            {
                sig_pairs[nodea].insert(nodeb);
                sig_pairs[nodeb].insert(nodea);
            } else{
                nonsig_pairs[nodea].insert( nodeb);
                nonsig_pairs[nodeb].insert( nodea );
            }
        }
    }

    ofstream outf( "node_sigP_stat.txt");
    outf<<"Node\tSigCorPair_n\tNonsigCorPair_n"<<endl;
    for ( map<string, string >::iterator ite = node_celltype.begin(); ite != node_celltype.end(); ++ite )
    {
        string node = ite->first;
        int sn = 0;
        int nn = 0;
        if ( sig_pairs.find(node) != sig_pairs.end() )
            sn += (int)sig_pairs[node].size();
            
        if ( nonsig_pairs.find(node) != nonsig_pairs.end() )
        {
            nn += (int)nonsig_pairs[node].size();
        }
        outf<<node<<"\t"<<sn<<"\t"<<nn<<endl;
    }
    outf.close();
}

void filter_node( map<string, set<string > > &sig_pairs,
    map<string, set<string > > &nonsig_pairs,
    map<string, string > &node_celltype,
    int sig_n_thr,
    map<string, set<string > > &sig_pairs_af,
    map<string, set<string > > &nonsig_pairs_af )
{
    set<string > outnode;
    set<string > innode;
    for ( map<string, set<string > >::iterator ite = sig_pairs.begin(); ite != sig_pairs.end(); ++ite )
    {
        string node = ite->first;
        int n = (int)ite->second.size();
        if ( n >= sig_n_thr )
        {
            innode.insert( node );
            sig_pairs_af[node] = ite->second;
        } else{
            outnode.insert( node );
        }
    }

    for ( map<string, set<string > >::iterator ite = nonsig_pairs.begin(); ite != nonsig_pairs.end(); ++ite )
    {
        string node = ite->first;
        if ( innode.find( node ) != innode.end() )
        {
            set<string > target_af;
            for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
            {
                if ( innode.find(*si ) != innode.end() )
                {
                    target_af.insert(*si);
                }
            }
            nonsig_pairs_af[node] = target_af;
        }
    }

    ofstream outf( "node_sigP_afterfilter_stat.txt");
    outf<<"Node\tSigCorPair_n\tNonsigCorPair_n"<<endl;
    for ( map<string, string >::iterator ite = node_celltype.begin(); ite != node_celltype.end(); ++ite )
    {
        string node = ite->first;
        int sn = 0;
        int nn = 0;
        
        if ( sig_pairs_af.find(node) != sig_pairs_af.end() )
            sn += (int)sig_pairs_af[node].size();
            
        if ( nonsig_pairs_af.find(node) != nonsig_pairs_af.end() )
        {
            nn += (int)nonsig_pairs_af[node].size();
        }
        if ( sn > 0 && nn > 0 )
            outf<<node<<"\t"<<sn<<"\t"<<nn<<endl;
    }
    outf.close();
}



void readin_ovlgene( string infile, map<string, map<string, set<string > > > &nodea_nodeb_ovlgene )
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
        vector<string > ps = parse_string( line, '\t' );
        string nodea = ps[0];
        string nodeb = ps[1];
        string ovlg_s = ps[7];
        vector<string > ovlg_ve = parse_string( ovlg_s, ',');
        set<string > ovlgs;
        for ( size_t i = 0; i < ovlg_ve.size(); ++i )
        {
            ovlgs.insert( ovlg_ve[i]);
        }
        nodea_nodeb_ovlgene[nodea][nodeb] = ovlgs;
        nodea_nodeb_ovlgene[nodeb][nodea] = ovlgs;
    }
    inf.close();

}

void readin_ligand_target2( string infile, map<string, map<string, map<string, set<string > > > > &nodeb_target_ligand_nodea )
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
        vector<string > ps = parse_string( line, '\t' );
        string nodea = ps[0];
        string nodeb = ps[1];
        string lt_s = ps[8];
        if ( lt_s == "NA" )
            continue;

        vector<string > lt_l1 = parse_string( lt_s, '|');
        
        for ( size_t i = 0; i < lt_l1.size(); ++i )
        {
            string lt_l1_s = lt_l1[i];
            vector<string > lt_l2 = parse_string( lt_l1_s, ':');
            string l = lt_l2[0];
            vector<string > lt_l3 = parse_string( lt_l2[1], ',');
            for ( size_t j = 0; j < lt_l3.size(); ++j )
            {
                string t = lt_l3[j];
                nodeb_target_ligand_nodea[nodeb][t][l].insert(nodea);
            }
        }
        
    }
    inf.close();
}

bool get_shared_ligand_nodea( map<string, set<string > > &ligand_nodea_A, map<string, set<string > > &ligand_nodea_B,
    map<string, set<string > > &shared_ligand_nodea, set<string > &shared_ligand )
{
    
    for ( map<string, set<string > >::iterator ite = ligand_nodea_A.begin(); ite != ligand_nodea_A.end(); ++ite )
    {
        string ligand = ite->first;
        if ( ligand_nodea_B.find( ligand ) == ligand_nodea_B.end() )
        {
            continue;
        } else
        {
            shared_ligand.insert( ligand );
        }

        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            if ( ligand_nodea_B[ligand].find( *si ) != ligand_nodea_B[ligand].end() )
            {
                shared_ligand_nodea[ligand].insert( *si );
               
            }
        }
    }
    if ( shared_ligand_nodea.empty() && shared_ligand.empty() )
        return false;
    else
        return true;
    
}

void integrate_ovlgene_ligandtarget( map<string, map<string, pair<double,double> > > &cor_matrix,
    map<string, map<string, set<string > > > &nodea_nodeb_ovlgene,
    map<string, map<string, map<string, set<string > > > > &nodeb_target_ligand_nodea,
    string outfile, string outfile2 ) 
{
    ofstream outf( outfile.data() );
    outf<<"nodea\tnodeb\tcor\tcorp\tSig\tn_ovlg\tn_suc_target\tligand_target_node\tn_suc_target_loose\tligand_target"<<endl;
    ofstream outf2( outfile2.data() );
    outf2<<"node\tSig\tn_agg_ovlg\tn_agg_suc_target\tn_agg_suc_target_loose\tratio_agg_suc_target\tratio_agg_suc_target_loose"<<endl;

    map<string, int > Sig_agg_ovlg;
    map<string, int > Sig_agg_suc_target;
    map<string, int > Sig_agg_suc_target_loose;
    map<string, int > NonSig_agg_ovlg;
    map<string, int > NonSig_agg_suc_target;
    map<string, int > NonSig_agg_suc_target_loose;

    for ( map<string, map<string, set<string > > >::iterator ite = nodea_nodeb_ovlgene.begin(); ite != nodea_nodeb_ovlgene.end(); ++ite )
    {
        string nodea = ite->first;
        int sig_agg_ovlg = 0;
        int sig_agg_suc_target = 0;
        int sig_agg_suc_target_loose = 0;
        int nonsig_agg_ovlg = 0;
        int nonsig_agg_suc_target = 0;
        int nonsig_agg_suc_target_loose = 0;

        for ( map<string, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string nodeb = si->first;

            pair<double, double > cors = cor_matrix[nodea][nodeb];

            string sig = "NonSig";
            if ( cors.first >= 0.5 & cors.second < 0.001 )
                sig = "Sig";

            int n_ovlg = (int)si->second.size();

            if ( *si->second.begin( ) == "NONE" )
            {
                outf<<nodea<<"\t"<<nodeb<<"\t"<<cors.first<<"\t"<<cors.second<<"\t"<<sig<<"\t"<<n_ovlg<<"\t0\tNA\t0\tNA"<<endl;
                continue;
            }

            if ( sig == "Sig" )
            {
                sig_agg_ovlg += n_ovlg;
            } else
                nonsig_agg_ovlg += n_ovlg;

            map<string, map<string, set<string > > > suc_target;
            map<string, set<string > > suc_target_loose;
            if ( (nodeb_target_ligand_nodea.find( nodea ) != nodeb_target_ligand_nodea.end()) && 
                (nodeb_target_ligand_nodea.find( nodeb ) != nodeb_target_ligand_nodea.end()) )
            {
                for ( set<string >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
                {
                    string t = *ti;
                    if ( nodeb_target_ligand_nodea[nodea].find( t ) != nodeb_target_ligand_nodea[nodea].end() &&
                        (nodeb_target_ligand_nodea[nodeb].find( t ) != nodeb_target_ligand_nodea[nodeb].end()) )
                    {
                        
                        map<string, set<string > > shared_ligand_nodea;
                        set<string > shared_ligand;
                        if ( get_shared_ligand_nodea(nodeb_target_ligand_nodea[nodea][t], nodeb_target_ligand_nodea[nodeb][t], shared_ligand_nodea, shared_ligand ) )
                        {
                            if ( !shared_ligand_nodea.empty() )
                                suc_target[t] = shared_ligand_nodea;
                            suc_target_loose[t] = shared_ligand;
                        } 
                    }
                    
                }
            }

            outf<<nodea<<"\t"<<nodeb<<"\t"<<cors.first<<"\t"<<cors.second<<"\t"<<sig<<"\t"<<n_ovlg;
            if ( !suc_target.empty() )
            {
                int n_suc_target = (int)suc_target.size();

                if ( sig == "Sig" )
                {
                    sig_agg_suc_target += n_suc_target;
                } else
                {
                    nonsig_agg_suc_target += n_suc_target;
                }

                map<string, set<string > > lt;
                map<string, set<string > > l_na;
                for ( map<string, map<string, set<string > > >::iterator ti = suc_target.begin(); ti != suc_target.end(); ++ti )
                {
                    string t = ti->first;
                    for ( map<string, set<string > >::iterator tti = ti->second.begin(); tti != ti->second.end(); ++tti )
                    {
                        string l = tti->first;
                        lt[l].insert(t);
                        for ( set<string >::iterator ni = tti->second.begin(); ni != tti->second.end(); ++ni )
                        {
                            string na = *ni;
                            l_na[l].insert(na);
                        }
                    }
                } 


                outf<<"\t"<<n_suc_target;
                for ( map<string, set<string > >::iterator lti = lt.begin(); lti != lt.end(); ++lti )
                {
                    string l = lti->first;
                    if ( lti == lt.begin() )
                    {
                        outf<<"\t";
                    } else{
                        outf<<"|";
                    }
                    outf<<lti->first;
                    for ( set<string>::iterator slti = lti->second.begin(); slti != lti->second.end(); ++slti )
                    {
                        if ( slti == lti->second.begin() )
                            outf<<":";
                        else
                            outf<<",";

                        outf<<*slti;
                    }
                    for ( set<string>::iterator slna = l_na[l].begin(); slna != l_na[l].end(); ++slna )
                    {
                        if ( slna == l_na[l].begin() )
                        {
                            outf<<":";
                        } else
                            outf<<",";
                        outf<<*slna;
                    }
                }
                
                
            } else {
                outf<<"\t0\tNA";
            }
            if ( !suc_target_loose.empty() )
            {
                int n_suc_target_loose = (int)suc_target_loose.size();

                if ( sig == "Sig" )
                {
                    sig_agg_suc_target_loose += n_suc_target_loose;
                } else
                {
                    nonsig_agg_suc_target_loose += n_suc_target_loose;
                }

                outf<<"\t"<<n_suc_target_loose;
                map<string, set<string > > lt;
                for ( map<string, set<string > >::iterator ti = suc_target_loose.begin(); ti != suc_target_loose.end(); ++ti )
                {
                    string t = ti->first;
                    for ( set<string >::iterator lti = ti->second.begin(); lti != ti->second.end(); ++lti )
                    {
                        lt[*lti].insert( t );
                    }
                }

                for ( map<string, set<string > >::iterator lti = lt.begin(); lti != lt.end(); ++lti )
                {
                    string l = lti->first;
                    if ( lti == lt.begin() )
                    {
                        outf<<"\t";
                    } else{
                        outf<<"|";
                    }
                    outf<<lti->first;
                    for ( set<string>::iterator slti = lti->second.begin(); slti != lti->second.end(); ++slti )
                    {
                        if ( slti == lti->second.begin() )
                            outf<<":";
                        else
                            outf<<",";

                        outf<<*slti;
                    }
                }

            } else{
                outf<<"\t0\tNA";
            }

            outf<<endl;

            
        }

        Sig_agg_ovlg[nodea] = sig_agg_ovlg;
        Sig_agg_suc_target[nodea] = sig_agg_suc_target;
        Sig_agg_suc_target_loose[nodea] = sig_agg_suc_target_loose;
        NonSig_agg_ovlg[nodea] = nonsig_agg_ovlg;
        NonSig_agg_suc_target[nodea] = nonsig_agg_suc_target;
        NonSig_agg_suc_target_loose[nodea] = nonsig_agg_suc_target_loose;


        double sig_agg_suc_target_r = 0;
        double sig_agg_suc_target_loose_r = 0;
        if ( sig_agg_ovlg > 0 )
        {
            sig_agg_suc_target_r = sig_agg_suc_target * 1.0/sig_agg_ovlg;
            sig_agg_suc_target_loose_r = sig_agg_suc_target_loose * 1.0 / sig_agg_ovlg;
        } 
        double nonsig_agg_suc_target_r = 0;
        double nonsig_agg_suc_target_loose_r = 0;
        if ( nonsig_agg_ovlg > 0 )
        {
            nonsig_agg_suc_target_r = nonsig_agg_suc_target * 1.0/nonsig_agg_ovlg;
            nonsig_agg_suc_target_loose_r = nonsig_agg_suc_target_loose * 1.0 / nonsig_agg_ovlg;
        } 

        outf2<<nodea<<"\tSig\t"<<sig_agg_ovlg<<"\t"<<sig_agg_suc_target<<"\t"<<sig_agg_suc_target_loose<<"\t"<<sig_agg_suc_target_r<<"\t"<<sig_agg_suc_target_loose_r<<endl;
        outf2<<nodea<<"\tNonSig\t"<<nonsig_agg_ovlg<<"\t"<<nonsig_agg_suc_target<<"\t"<<nonsig_agg_suc_target_loose<<"\t"<<nonsig_agg_suc_target_r<<"\t"<<nonsig_agg_suc_target_loose_r<<endl;
        
    }
}

void readin_ligand_target( string infile, 
    map<string, map<string, pair<double, double > > > &nodea_nodeb_ra_rb,
    map<string, map<string, pair<double, double > > > &nodeb_nodea_ra_rb )
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
        vector<string > ps = parse_string( line, '\t' );
        string nodea = ps[0];
        string nodeb = ps[1];
        double ra = atof( ps[6].c_str() );
        double rb = atof( ps[7].c_str() );

        nodea_nodeb_ra_rb[nodea][nodeb] = make_pair( ra, rb );
        nodeb_nodea_ra_rb[nodeb][nodea] = make_pair( ra, rb );
        
    }
    inf.close();

}

double mean_r( vector<double > &ve)
{
    if ( ve.empty() )
        return 0;
    double a = 0;
    int s = (int)ve.size();
    for ( size_t i = 0; i < ve.size(); ++i )
        a += ve[i];
    double m = a / s;
    return m;
}

void integrate( map<string, set<string > > &sig_pairs_af,
    map<string, set<string > > &nonsig_pairs_af,
    map<string, map<string, pair<double, double > > > &nodea_nodeb_ra_rb,
    map<string, map<string, pair<double, double > > > &nodeb_nodea_ra_rb,
    string outfile1, string outfile2 )
{
    ofstream outf1( outfile1.data() );

    outf1<<"TargetNode\tSigPair_n\tNonsigPair_n\tSigPair_target_meanratio\tNonsigPair_target_meanratio"<<endl;
    
    for ( map<string, map<string, pair<double, double > > >::iterator ite = nodeb_nodea_ra_rb.begin();
        ite != nodeb_nodea_ra_rb.end(); ++ite )
    {
        string nodeb = ite->first;
        if ( sig_pairs_af.find( nodeb ) == sig_pairs_af.end() )
            continue;
            
        vector<double> sig_target_r_ve;
        vector<double> nonsig_target_r_ve;
    //    int sig_lig_n = 0;
    //    int nonsig_lig_n = 0;
        for ( set<string >::iterator si = sig_pairs_af[nodeb].begin(); si != sig_pairs_af[nodeb].end(); ++si )
        {
            string nodea = *si;
            if ( nodeb_nodea_ra_rb[nodeb].find( nodea ) == nodeb_nodea_ra_rb[nodeb].end() )
            {
                cout<<"warning cannot find nodea "<<nodea<<" for nodeb "<<nodeb<<endl; exit(1);
            }
        /*    double ra = nodeb_nodea_ra_rb[nodeb][nodea].first;
            if ( ra > 0 )
            {
                sig_lig_n += 1;
            } */
            double rb = nodea_nodeb_ra_rb[nodea][nodeb].second;
           
            sig_target_r_ve.push_back(rb);
        }
        for ( set<string >::iterator si = nonsig_pairs_af[nodeb].begin(); si != nonsig_pairs_af[nodeb].end(); ++si )
        {
            string nodea = *si;
            if ( nodeb_nodea_ra_rb[nodeb].find( nodea ) == nodeb_nodea_ra_rb[nodeb].end() )
            {
                cout<<"warning cannot find nodea "<<nodea<<" for nodeb "<<nodeb<<endl; exit(1);
            }
         /*   double ra = nodeb_nodea_ra_rb[nodeb][nodea].first;
            if ( ra > 0 )
            {
                nonsig_lig_n += 1;
            }  */
            double rb = nodea_nodeb_ra_rb[nodea][nodeb].second;
            
            nonsig_target_r_ve.push_back(rb);
        }

        int SigPair_n = (int)sig_pairs_af[nodeb].size();
        int NonsigPair_n = (int)nonsig_pairs_af[nodeb].size();
     //   double Sig_ratio = (sig_lig_n *1.0 )/SigPair_n;
     //   double Nonsig_ratio = (nonsig_lig_n*1.0)/NonsigPair_n;
        double Sig_target_meanratio = mean_r( sig_target_r_ve );
        double Nonsig_target_meanratio = mean_r( nonsig_target_r_ve );

        outf1<<nodeb<<"\t"<<SigPair_n<<"\t"<<NonsigPair_n<<"\t"<< Sig_target_meanratio<<"\t"<<Nonsig_target_meanratio<<endl;
    }
    outf1.close();
    
    ofstream outf2( outfile2.data() ); 

    outf2<<"SenderNode\tSigPair_n\tNonsigPair_n\tSig_target_meanratio\tNonsig_target_meanratio"<<endl;

    for ( map<string, map<string, pair<double, double > > >::iterator ite = nodea_nodeb_ra_rb.begin();
        ite != nodea_nodeb_ra_rb.end(); ++ite )
    {
        string nodea = ite->first;
        if ( sig_pairs_af.find( nodea ) == sig_pairs_af.end() )
            continue;
        vector<double> sig_target_r_ve;
        vector<double> nonsig_target_r_ve;
        
        for ( set<string >::iterator si = sig_pairs_af[nodea].begin(); si != sig_pairs_af[nodea].end(); ++si )
        {
            string nodeb = *si;
            if ( nodea_nodeb_ra_rb[nodea].find( nodeb ) == nodea_nodeb_ra_rb[nodea].end() )
            {
                cout<<"warning cannot find nodeb "<<nodeb<<" for nodea "<<nodea<<endl; exit(1);
            }
            double rb = nodea_nodeb_ra_rb[nodea][nodeb].second;
            sig_target_r_ve.push_back(rb);
        }

        for ( set<string >::iterator si = nonsig_pairs_af[nodea].begin(); si != nonsig_pairs_af[nodea].end(); ++si )
        {
            string nodeb = *si;
            if ( nodea_nodeb_ra_rb[nodea].find( nodeb ) == nodea_nodeb_ra_rb[nodea].end() )
            {
                cout<<"warning cannot find nodeb "<<nodeb<<" for nodea "<<nodea<<endl; exit(1);
            }
            double rb = nodea_nodeb_ra_rb[nodea][nodeb].second;
            nonsig_target_r_ve.push_back(rb);
        }
        int SigPair_n = (int)sig_pairs_af[nodea].size();
        int NonsigPair_n = (int)nonsig_pairs_af[nodea].size();
        double Sig_target_meanratio = mean_r( sig_target_r_ve );
        double Nonsig_target_meanratio = mean_r( nonsig_target_r_ve );

        outf2<<nodea<<"\t"<<SigPair_n<<"\t"<<NonsigPair_n<<"\t"<<Sig_target_meanratio<<"\t"<<Nonsig_target_meanratio<<endl;

    }
    outf2.close();  
}



int main(int argc, char* argv[] )
{
    if ( argc == 1 )
    {
        cout<<"combine link cor and metrics "<<endl;
    //    cout<<"Usage: prog linkfile vertfile mode inltfile| outfile"<<endl;
        cout<<"Usage: prog linkfile ovlgfile ltfile outfile outfile2"<<endl;
        exit(1);
    }
/*    string infile1 = argv[1];
    string infile2 = argv[2];
    int mode = atoi( argv[3]);

    double pthr = 0.001;
    double corthr = 0.5;

    map<string, map<string, pair<double,double> > > cor_matrix;
    readinlinkcor( infile1, cor_matrix );

    map<string, string > node_celltype;
    readinvert( infile2, node_celltype );


    map<string, set<string > > sig_pairs;
    map<string, set<string > > nonsig_pairs;
    seperate( cor_matrix,  node_celltype, sig_pairs, nonsig_pairs, pthr, corthr );

    map<string, set<string > > sig_pairs_af;
    map<string, set<string > > nonsig_pairs_af;
    int sig_n_thr = 5;
    filter_node( sig_pairs, nonsig_pairs, node_celltype, sig_n_thr, sig_pairs_af, nonsig_pairs_af );

    if ( mode == 1 )
    {
        string inltfile = argv[4];
        string outfile1 = argv[5];
        string outfile2 = argv[6];

        map<string, map<string, pair<double, double > > > nodea_nodeb_ra_rb;
        map<string, map<string, pair<double, double > > > nodeb_nodea_ra_rb;
        readin_ligand_target( inltfile, nodea_nodeb_ra_rb, nodeb_nodea_ra_rb );

        integrate( sig_pairs_af, nonsig_pairs_af, nodea_nodeb_ra_rb, nodeb_nodea_ra_rb, outfile1, outfile2 );
    }  */

    string inlinkfile = argv[1];
    string inovlgfile = argv[2];
    string inltfile = argv[3];
    string outfile = argv[4];
    string outfile2 = argv[5];

    map<string, map<string, pair<double,double> > > cor_matrix;
    readinlinkcor( inlinkfile, cor_matrix );

    map<string, map<string, set<string > > > nodea_nodeb_ovlgene;
    readin_ovlgene( inovlgfile, nodea_nodeb_ovlgene );

    map<string, map<string, map<string, set<string > > > > nodeb_target_ligand_nodea;
    readin_ligand_target2( inltfile, nodeb_target_ligand_nodea );

    cout<<"Integrate"<<endl;

    integrate_ovlgene_ligandtarget( cor_matrix, nodea_nodeb_ovlgene, nodeb_target_ligand_nodea, outfile, outfile2 );

    return 1;

}


