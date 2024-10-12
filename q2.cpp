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

vector<vector<double>> getProfile(vector<string> &motifs, int exceptId){
    vector<vector<double>> profile(motifs[0].size(),vector<double>(4,1.0));
    for(int i=0;i<motifs.size();i++){
        if(i == exceptId){
            continue;
        }
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

int getRandomIndex(const std::vector<double>& weights) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    int randomIndex = dist(gen);

    return randomIndex;
}

void getWeightedProbabilityUpdatedMotif(vector<string> &motifs,vector<vector<double>> &profile,string &consensus,vector<string> &dna,int motifId){

    string chosenMotif = dna[motifId];

    vector<double> probs;

    int mlen = dna[motifId].size()-consensus.size();

    for(int x=0;x<=mlen;x++){
        double mscore = getProbabilityMotifScore(profile,dna[motifId],x);

        probs.push_back(mscore);
    }
    int ind = getRandomIndex(probs);
    motifs[motifId] = dna[motifId].substr(ind,consensus.size());        
}


struct MotifResult{
    vector<string> motifs;
    int score;
};



MotifResult GibbsSampling(vector<string> &dna,int k, int t,int r){
    
    int n = dna[0].size();
    MotifResult bestMotifResult;
    bestMotifResult.score = INT_MAX;

    auto motifs = getRandomMotif(n,k,dna);
    
    for(int i=0;i<r;i++){
        int chosenMotifInd = rand()%t;
        auto profile = getProfile(motifs,chosenMotifInd);
        string consensus = getConsensus(profile);
        getWeightedProbabilityUpdatedMotif(motifs,profile,consensus,dna,chosenMotifInd);
        
        int motifScore = getMotifScore(motifs,consensus);

        if(motifScore<bestMotifResult.score){
            bestMotifResult.score = motifScore;
            bestMotifResult.motifs = motifs;
        }
    }

    return bestMotifResult;
}

void hw2q2(string filename){

    auto stripped = split(filename,'.')[0];
    auto id = split(stripped,'_').back();
    int fnum = atoi(id.c_str());

    string s = read_string_from_file(filename);
    vector<string> sp = split(s,'\n');
    string nums = sp[0];
    auto splitnums = split(nums,' ');
    int k = atoi(splitnums[0].c_str());
    int t = atoi(splitnums[1].c_str());
    int r = atoi(splitnums[2].c_str());
    vector<string> dna(sp.begin()+1,sp.end());

    MotifResult final_result;
    final_result.score = INT_MAX;
    for(int i=0;i<30;i++){
        
        
        auto res = GibbsSampling(dna,k,t,r);
        if(res.score<final_result.score){
            final_result.score = res.score;
            final_result.motifs = res.motifs;
        }
    }
    

    string output = "";

    for(auto x: final_result.motifs){
        output += x;
        output += '\n';
        cout<<x<<endl;
    }
    
    cout<<final_result.score;

    if (filename.find(string("test_")) != std::string::npos){
        cout<<"\nsaved\n";
        output_to_file(output,"sol_q2_t"+to_string(fnum)+".txt");
    }

}


void hw2q3(){

    vector<int> iterations = {1000, 2000, 10000};

    for(auto iter:iterations){
        string s = read_string_from_file("motif_dataset.txt");
        vector<string> sp = split(s,'\n');
        vector<string> dna(sp.begin(),sp.end());

        MotifResult final_result;
        final_result.score = INT_MAX;
        for(int i=0;i<30;i++){

            auto res = GibbsSampling(dna,15,10,iter);
            if(res.score<final_result.score){
                final_result.score = res.score;
                final_result.motifs = res.motifs;
            }
        }

        auto profile = getProfile(final_result.motifs,-1);
        string consensus = getConsensus(profile);
        
        cout<<consensus<<endl;
        
        cout<<final_result.score<<endl<<endl;
    }

}

#include <filesystem>
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    srand(time(0));
    
    if(argc==2 ){
        string cmd = argv[1];
        string file = string(argv[1]);
        hw2q2(file);
    }  
    else{
        hw2q3();
    }  
    return 0;
}