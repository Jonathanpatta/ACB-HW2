#include <bits/stdc++.h>
using namespace std;

void output_to_file(string input,string filename){
    std::ofstream out(filename);
    out << input;
    out.close();
}

std::string read_string_from_file(const std::string &file_path) {

    const std::ifstream input_stream(file_path);
    if (input_stream.fail()) {
        throw std::runtime_error("Failed to open file");
    }
    string s;
    if(input_stream){
        std::ostringstream buffer;
        buffer << input_stream.rdbuf();
        s = buffer.str();
    }
    return s;
}

vector<string> split(string s,char c){

    std::stringstream test(s);
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(test, segment, c))
    {
       seglist.push_back(segment);
    }
    return seglist;
}

vector<string> getRandomMotif(int n, int k, vector<string> &dna){
    vector<string> motifs;
    for(auto x: dna){
        int ind = rand()%(n-(k+1));
        motifs.push_back(x.substr(ind,k));
    }
    return motifs;
}

vector<vector<double>> getProfile(vector<string> &motifs){
    vector<vector<double>> profile(motifs[0].size(),vector<double>(4,1.0));
    for(int i=0;i<motifs.size();i++){
        for(int j=0;j<motifs[i].size();j++){
            switch (motifs[i][j])
            {
            case 'A':
                profile[j][0]+=1;
                break;
                
            case 'C':
                profile[j][1]+=1;
                break;
                
            case 'G':
                profile[j][2]+=1;
                break;
                
            case 'T':
                profile[j][3]+=1;
                break;
            
            default:
                break;
            }
        }
    }
    for(int i=0;i<profile.size();i++){
        profile[i][0]/=(double)(motifs.size()+4);
        profile[i][1]/=(double)(motifs.size()+4);
        profile[i][2]/=(double)(motifs.size()+4);
        profile[i][3]/=(double)(motifs.size()+4);
    }
    return profile;
}

string getConsensus(vector<vector<double>> &profile){
    string ans = "";

    for(auto x: profile){
        int index = std::distance(x.begin(), max_element(x.begin(),x.end()));
        if(index==0){
            ans += 'A';
        }
        if(index==1){
            ans += 'C';
        }
        if(index==2){
            ans += 'G';
        }
        if(index==3){
            ans += 'T';
        }
    }

    return ans;

}

double getProbabilityMotifScore(vector<vector<double>> &profile,string &s,int start){
    double ans = 1;
    for(int i=0;i<profile.size();i++){
        switch (s[i+start])
        {
        case 'A':
            ans*=profile[i][0];
            break;
            
        case 'C':
            ans*=profile[i][1];
            break;
            
        case 'G':
            ans*=profile[i][2];
            break;
            
        case 'T':
            ans*=profile[i][3];
            break;
        
        default:
            break;
        }
    }
    return ans;
}

int getMotifScore(vector<string> &motifs, string &consensus){
    int score = 0;
    for(auto motif:motifs){
        for(int i=0;i<motif.size();i++){
            if(motif[i]!=consensus[i]) score++;
        }
    }
    return score;
}

void getProfileMostProbableMotif(vector<string> &motifs,vector<vector<double>> &profile,string &consensus,vector<string> &dna){

    for(int i=0;i<dna.size();i++){
        double score = 0.0;
        int start = -1;
        for(int x=0;x<=dna[i].size()-consensus.size();x++){
            double mscore = getProbabilityMotifScore(profile,dna[i],x);
            if(mscore>score){
                score = mscore;
                start = x;
            }
        }
        motifs[i] = dna[i].substr(start,consensus.size());
    }
}

struct MotifResult{
    vector<string> motifs;
    int score;
};

MotifResult RandomizedMotifSearch(vector<string> &dna,int k, int t){
    
    int n = dna[0].size();
    auto motifs = getRandomMotif(n,k,dna);
    MotifResult bestMotifResult;
    bestMotifResult.score = INT_MAX;
    
    while(true){
        auto profile = getProfile(motifs);
        string consensus = getConsensus(profile);

        getProfileMostProbableMotif(motifs,profile,consensus,dna);
        
        int motifScore = getMotifScore(motifs,consensus);

        if(motifScore<bestMotifResult.score){
            bestMotifResult.score = motifScore;
            bestMotifResult.motifs = motifs;
        }
        else{
            break;
        }
    }
    return bestMotifResult;
}

void hw2q1(string filename){

    auto stripped = split(filename,'.')[0];
    auto id = split(stripped,'_').back();
    int fnum = atoi(id.c_str());

    string s = read_string_from_file(filename);
    vector<string> sp = split(s,'\n');
    string nums = sp[0];
    auto splitnums = split(nums,' ');
    int k = atoi(splitnums[0].c_str());
    int t = atoi(splitnums[1].c_str());
    vector<string> dna(sp.begin()+1,sp.end());


    MotifResult r;
    r.score = INT_MAX;
    for(int i=0;i<1500;i++){        
        auto res = RandomizedMotifSearch(dna,k,t);
        if(res.score<r.score){
            r.score = res.score;
            r.motifs = res.motifs;
        }
    }
    string output = "";

    for(auto x: r.motifs){
        output += x;
        output += '\n';
        cout<<x<<endl;
    }
    
    cout<<r.score;
    if (filename.find(string("test_")) != std::string::npos){
        cout<<"\nsaved\n";
        output_to_file(output,"sol_q1_t"+to_string(fnum)+".txt");
    }

}


void hw2q3(){

    vector<int> iterations = {1000, 10000, 100000};

    for(auto iter:iterations){
        string s = read_string_from_file("motif_dataset.txt");
        vector<string> sp = split(s,'\n');
        vector<string> dna(sp.begin(),sp.end());


        MotifResult r;
        r.score = INT_MAX;
        for(int i=0;i<iter;i++){        
            auto res = RandomizedMotifSearch(dna,15,10);
            if(res.score<r.score){
                r.score = res.score;
                r.motifs = res.motifs;
            }
        }
        // string output = "";

        // for(auto x: r.motifs){
        //     output += x;
        //     output += '\n';
        //     cout<<x<<endl;
        // }
        
        auto profile = getProfile(r.motifs);
        string consensus = getConsensus(profile);
        
        cout<<consensus<<endl;
        cout<<r.score<<endl<<endl;
    }

}



#include <filesystem>
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    srand(time(0));
    
    if(argc==2 ){
        string cmd = argv[1];
        string file = string(argv[1]);
        hw2q1(file);
    } 
    else{
        hw2q3();
    }   
    return 0;
}