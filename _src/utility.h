#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <ctime>
#include <string>
#include <vector>

using namespace std;

// A loose collection of utility functions for inclusion to main

void print_banner(const char* fname){
    ifstream ifs(fname, ifstream::in);

    char c = ifs.get();
    while (ifs.good()) {
        std::cout << c;
        c = ifs.get();
    }
    ifs.close();
}

void try_mkdir(const std::string& fname) {
    if (filesystem::create_directory(std::string(fname))) {
        if (errno != EEXIST)
            cerr<<"mkdir error attempting to create "<< fname << ":" << errno;
    }
}

/// Copies .mol file (e.g. to log directory).
void save_mol_file(const std::string& path_to_file,const std::string& out_path) {
    
    cout << "[ Input ] Storing input in directory "<<out_path<<"..."<<endl;
    std::filesystem::path infile_path = string(path_to_file);
    std::filesystem::path outfile_path = out_path; // infile_path.filename() Returns "MoleculeName.mol"

    try{
        std::filesystem::copy_file(infile_path,outfile_path);
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
}
/// Moves .mol file (e.g. from log to output directory).
string move_mol_file(const string& path_to_file,const string& out_dir, const string& tag) {
    cout << "[ Input ] Copying " << std::filesystem::path(path_to_file).filename() << " to directory "<<out_dir<<"..."<<endl;
    std::filesystem::path infile_path = string(path_to_file);
    std::filesystem::path outfile_path = out_dir + tag + ".mol";
    try{
        std::filesystem::copy_file(infile_path,outfile_path);
        std::filesystem::remove(infile_path);
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
    try{
        std::filesystem::remove(infile_path);
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }    
    return outfile_path;
}
string move_log_file(const string& path_to_file,const string& out_dir, const string& tag) {
    cout << "[ Log ] Copying " << std::filesystem::path(path_to_file).filename() << " to directory "<<out_dir<<"..."<<endl;
    std::filesystem::path infile_path = string(path_to_file);
    std::filesystem::path outfile_path = out_dir + tag + ".log";
    try{
        std::filesystem::copy_file(infile_path,outfile_path);
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }    
    try{
        std::filesystem::remove(infile_path);
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
    return outfile_path;
}

bool string_in_file(vector<string> FileContent, string outdir){
    for (size_t n = 0; n < FileContent.size(); n++){
        if (FileContent[n] == outdir){
            return true;
        }
    }    
    return false;
}

int reserve_output_name(string &outdir,string &tag){
    ifstream f_read;
    string fname = "output/reserved_output_names.txt";
    f_read.open(fname);    
    vector<string> FileContent;
	string comment = "//";
	string curr_key;
    // Store file content
	while (!f_read.eof() && f_read.is_open()) {
		string line;
		getline(f_read, line);
		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 2, "")) continue;
		FileContent.push_back(line);
	}

    int count = 0;
    outdir = "";
    while(outdir == ""||string_in_file(FileContent,outdir)){
            count++;
            outdir = "output/__Molecular/"+tag+"_"+to_string(count)+"/";
            if(count > 999){
                cerr << "Directory naming loop failsafe triggered in get_file_names() of main.cpp" << endl;
                return 1;
            }
        }
    f_read.close();
    // Reserve the name in the file
    cout << "Reserving folder name: '" << fname << "'" << endl;
    ofstream f_write;
    f_write.open(fname,fstream::app); 
    f_write << endl <<outdir;        
    // Ensure we can make the directory
    try_mkdir(outdir);
    std::filesystem::remove(outdir);
    f_write.close();
    return 0;
}

int get_file_names(string &infile_, string &tag, string &tmp_logfile, string &tmp_molfile, string&outdir) {
    // Takes infile of the form "DIR/Lysozyme.mol"
    // Stores "Lysozyme" in tag, "output/log/run_Lysozyme" in logfile, and "output/log/Lysozyme_[timestamp].mol" in tmp_molfile
    
    if (!std::filesystem::exists(infile_)){
        // Allow for user to not include extension.
        infile_ += ".mol";
        if (!std::filesystem::exists(infile_)){
            cerr <<"File "<< infile_ <<" not found."<<endl;
            return 1;
        }
    }

    size_t tagstart = infile_.rfind('/');
    size_t tagend = infile_.rfind('.');
    tagstart = (tagstart==string::npos) ? 0 : tagstart + 1;// Exclude leading slash
    tagend = (tagend==string::npos) ? infile_.size() : tagend;
    tag = infile_.substr(tagstart, tagend-tagstart);

    time_t time_temp = std::time(nullptr); 
    tm time = *localtime(&time_temp);
    ostringstream time_tag;
    time_tag << std::put_time(&time, "%d-%m-%Y %H-%M-%S");
    
    // guarantee the existence of a folder structure
    try_mkdir("output");
    try_mkdir("output/log");
    try_mkdir("output/__Molecular");
    tmp_logfile = "output/log/run_" + tag + "_" + time_tag.str() + ".log";
    tmp_molfile = "output/log/mol_" + tag + "_" + time_tag.str() + ".mol"; 
    if (reserve_output_name(outdir,tag) == 1){return 1;}
    // check correct format
    string extension = infile_.substr(tagend);
    if (extension != ".mol") {
        cerr<<"This file is for coupled calculations. Please provide a .mol file similar to Lysozyme.mol"<<endl;
        return 1;
    }
    return 0;
}

struct CmdParser{
    CmdParser(int argc, const char *argv[]) {
        if (argc < 2) {
            std::cout << "Usage: solver path/to/molecular/in.mol [-sxhwc]" << std::endl;
            valid_input = false;
        }
        
        for (int a=2; a<argc; a++) {
            if (argv[a][0] != '-')
                continue;

            int i=1;
            while (argv[a][i]!='\0') {
                switch (argv[a][i]) {
                    case 's':
                        // recalculate
                        cout<<"\033[1;35mWarning:\033[0m Searching output folder for precalculated rates."<<endl;
                        cout<<"If these were calculated for a different X-ray energy, the results will be wrong!"<<endl;
                        recalc = false;
                        break;
                    case 'x':
                        // X - sections only
                        cout<<"\033[1;31mSolving for cross-sections only.\033[0m"<<endl;
                        solve_rate_eq = false;
                        break;
                    case 'h':
                        // Usage help.
                        cout<<"This is physics code, were you really expecting documentation?"<<endl;
                        cout<<"  -s Look for stored precalculated rate coefficients"<<endl;
                        cout<<"  -x Skip rate-equaton solving"<<endl;
                        break;
                    case 'w':
                        // Warranty.
                        cout<<"This program is provided as-is, with no warranty, explicit or implied."<<endl;
                        cout<<"It is a simplified model of XFEL plasma dynamics, however it should not"<<endl;
                        cout<<"be viewed as accurate under all conditions."<<endl;
                        exit(0);
                        break;
                    case 'c':
                        // Copyright.
                        print_banner("LICENSE");
                        exit(0);
                        break;

                    default:
                        cout<<"Flag '"<<argv[a][i]<<"' is not a recognised flag."<<endl;
                }
                i++;
            }
        }
    }
    bool recalc = true;
    bool valid_input = true;
    bool solve_rate_eq = true;
    bool load_data = false;
};