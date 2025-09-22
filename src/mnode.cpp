
// mnode.cpp

//BLock seeking by smOOth Moving
//BLOOM

// x distance, y distance
// inFile -> chr start end value
// outFile -> chr start end 

// max yi -> xi adjacent, yi adjacent yj -> xj adjacent, yj adjacet -> ... 
// reachable xi -> xj
// stop if (max y remained <= y distance)
// merge remained

// for each iteration , merge range
// 10 -- 30; 20 -- 50;  
// old merged : start -- end
// new merged : start -- end
// old and new overlap -> merge

// last merged range not output
// last merged range for each window; lastStart -- lastEnd -> add to new window
// last merged range == window --> new window == old window + 100

// default: xdis = 50, ydis = 0.05, win = 100, stepsize = 100  

// smooth, steady(high,low);  extend 50 bp
// link
//--------
// not merge <- search both side(peak); one side(expanding long)  

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <set>
#include <list>
#include <vector>
#include <cstring>
#include <cmath>
#include <dirent.h>
#include <sys/stat.h>


using namespace std;

typedef struct winval{
    int start;
    int pend;
    float value;
    int nnode;
/*    
    bool operator < (const struct winval &xval) const{

    }
*/

} WinDiv;

// peak, slop, step, \\  bottom
typedef struct mrange{
    int posStart;
    int posEnd;
    /*
    bool operator < (const struct mrange &xpos) const{
        if(posStart < xpos.posStart){
            return true;
        }
        if(posStart > xpos.posStart){
            return false;
        }
        
        if(posEnd < xpos.posEnd){
            return true;
        }
        return false;
    }
    */
} MnodeRan;

bool posSort(MnodeRan &a,MnodeRan &b){
    if(a.posStart < b.posStart){
        return true;
    }
    if(a.posStart > b.posStart){
        return false;
    }
    
    if(a.posEnd < b.posEnd){
        return true;
    }
    return false;
}


/*
double getValue(vector<WinDiv> &winVec,int winStart,int winEnd){
    double sumLen = 0.0,sumVal = 0.0;
    double kx2 = 0.0,sk = 0.0,kx = 0.0;
    for(int v = winStart; v <= winEnd; ++v){
        double tlen = winVec[v].pend - winVec[v].start + 1;
        sumLen += tlen;
        sumVal += winVec[v].value * tlen;
        double tk = tlen * winVec[v].value;
        kx2 += tk * winVec[v].value;
        kx += tk;
    }
    
    double wValue = sumVal / sumLen;
    return wValue;
}
*/

float wMean(vector<WinDiv> &winVec,int winStart,int winEnd,int wmet){
    double sumLen = 0.0,sumVal = 0.0;
    
    if(wmet == 1){
        for(int v = winStart; v <= winEnd; ++v){
            double tlen = winVec[v].pend - winVec[v].start + 1;
            sumLen += tlen;
            sumVal += winVec[v].value * tlen;
        }
    }else{
        for(int v = winStart; v <= winEnd; ++v){
            sumLen += winVec[v].nnode;
            sumVal += winVec[v].value * winVec[v].nnode;
        }
    }
    float wValue = sumVal / sumLen;
    return wValue;
}

static int saveCount = 0;
void storeValue(vector<WinDiv> &winVec,string &ass,string &tchr,int winStart,int winEnd,int wmet,ofstream &ofh){
    double sumLen = 0.0,sumVal = 0.0;
    //double kx2 = 0.0,kx = 0.0;
    int tlNode = 0;
    if(winStart == winEnd){
        ofh<<ass<<"\t"<<tchr<<"\t"<<winVec[winStart].start<<"\t"<<winVec[winEnd].pend<<"\t"<<fixed<<setprecision(5)<<winVec[winStart].value<<"\t"<<winVec[winStart].nnode<<endl;
    }else{
        float wValue = wMean(winVec,winStart,winEnd,wmet);
        for(int v = winStart; v <= winEnd; ++v){
            tlNode += winVec[v].nnode;
        }
        ofh<<ass<<"\t"<<tchr<<"\t"<<winVec[winStart].start<<"\t"<<winVec[winEnd].pend<<"\t"<<fixed<<setprecision(5)<<wValue<<"\t"<<tlNode<<endl; 
    }
    ++saveCount;
}

// max(dy,baseValue)
void mergeNode(vector<WinDiv> &winVec,string &ass,string &tchr,int dx,float dy,float baseValue,int plateau,int wmet,ofstream &ofh,int &lpos,int &rpos){
    vector<char> flag(winVec.size(),'0');
    vector<char> maxFlag(winVec.size(),'0');
    int oneCount = 0;
    int n = winVec.size();
    
    float maxValue = 0.0f;
    int k = 0;
    vector<MnodeRan> mran;
    //unordered_set<int> rangeSet;
    
    while(oneCount < n){
        
        bool fir = false;
        for(int i = 0; i < n; ++i){
            if(maxFlag[i] == '0'){
                if(! fir){
                    k = i;
                    maxValue = winVec[i].value;
                }else{
                    if(winVec[i].value > maxValue){
                        k = i;
                        maxValue = winVec[i].value;
                    }
                }
                fir = true;
            }
        }
        //
        if(maxValue > baseValue){
            maxFlag[k] = '1';
            ++oneCount;
            float upLimit,lowLimit;
            
            upLimit = maxValue + dy;
            lowLimit = maxValue - dy;
            
            int xend = winVec[k].pend + dx;
            int xstart = winVec[k].start - dx;
            int xlen = winVec[k].pend - winVec[k].start + 1;
            int rbound = k,lbound = k;
            bool rflag = false, lflag = false;
            //bool rmark = false, lmark = false;
            //int rmark = 0, lmark = 0;
            //float rlen = 0.0f;
            //float llen = 0.0f;
            // right
            int rpre = k,lpre = k;
            float rpVal = 0.0f,lpVal = 0.0f;
            for(int j = k + 1; j < n; ++j){
                if(winVec[j].start <= xend){
                    
                    
                    if(winVec[j].value >= lowLimit && winVec[j].value <= upLimit){
                        xend = winVec[j].pend + dx;                        
                        if(flag[j] == '0'){
                            rpre = rbound;
                            //
                            //maxFlag[j] = '1';
                            //++oneCount;
                            //
                            if(rpre > k){
                                float tval = wMean(winVec,rpre + 1,j,wmet);
                                float diff = tval - rpVal;
                                if(diff < 0){
                                    diff = -diff;
                                }
                                if(diff <= dy){
                                    rbound = j;
                                    //flag[j] = '1';
                                    //
                                    rflag = true;
                                }else{
                                    break;
                                }
                            }else{
                                if(winVec[j].start - winVec[rpre].pend > plateau){
                                    if(xlen >= plateau){
                                        // contain gap or another block
                                        break;
                                    }else{
                                        int ntLen = winVec[j].pend - winVec[j].start + 1;
                                        if(ntLen >= plateau){
                                            break;
                                        }else{
                                            if(j - rpre > 1){
                                                int mid = (winVec[j].start + winVec[rpre].pend) / 2;
                                                int rnum = 0,lnum = 0;
                                                int rStart = 0;
                                                for(int u = rpre + 1; u < j; ++u){
                                                    if(winVec[u].pend <= mid){
                                                        ++lnum;
                                                    }else if(winVec[u].start >= mid){
                                                        if(rnum == 0){
                                                            rStart = u;
                                                        }
                                                        ++rnum;
                                                    }
                                                }
                                                if(lnum == 0){
                                                    if(rnum == 0){
                                                        // cross middle
                                                        rbound = j;
                                                        //flag[j] = '1';
                                                        //
                                                        rpVal = wMean(winVec,rpre,j,wmet);
                                                        rflag = true;
                                                    }else{
                                                        break;
                                                    }
                                                }else{
                                                    break;
                                                    /*
                                                    if(rnum == 0){
                                                        break;
                                                    }else{
                                                        float lw = wMean(winVec,rpre + 1,rpre + lnum);
                                                        float rw = wMean(winVec,rStart,j - 1);
                                                        float mdiff = rw - lw;
                                                        if(mdiff < 0){
                                                            mdiff = -mdiff;
                                                        }
                                                        if(mdiff <= dy){
                                                            rbound = j;
                                                            //flag[j] = '1';
                                                            //
                                                            rpVal = wMean(winVec,rpre,j);
                                                            rflag = true;
                                                        }
                                                    }
                                                    */
                                                }
                                            }else{
                                                // gap
                                                break;
                                            }
                                        }
                                    }
                                }else{
                                    float tval = wMean(winVec,rpre,j,wmet);
                                    if(xlen >= plateau){
                                        float diff = winVec[rpre].value - tval;
                                        if(diff < 0){
                                            diff = -diff;
                                        }
                                        if(diff <= dy){
                                            rbound = j;
                                            //flag[j] = '1';
                                            rpVal = tval;
                                            rflag = true;
                                        }else{
                                            break;
                                        }
                                    }else{
                                        int ntLen = winVec[j].pend - winVec[j].start + 1;
                                        if(ntLen >= plateau){
                                            float diff = winVec[j].value - tval;
                                            if(diff < 0){
                                                diff = -diff;
                                            }
                                            if(diff <= dy){
                                                rbound = j;
                                                //flag[j] = '1';
                                                rpVal = tval;
                                                rflag = true;
                                            }else{
                                                break;
                                            }
                                        }else{
                                            rbound = j;
                                            //flag[j] = '1';
                                            rpVal = tval;
                                            rflag = true;
                                        }
                                    }
                                    
                                }
                            }
                            
                        }else{
                            break;
                        }
                    }else{
                        if(flag[j] == '1'){
                            break;
                        }
                    }
                }else{
                    break;
                }
            }
/*            
            if(rflag){
                if(rbound - rpre > 1){
                    if(winVec[rbound].start - winVec[rpre].pend > plateau){
                        if(winVec[rbound].pend - winVec[rbound].start < plateau){
                            flag[rbound] = '0';
                        }
                        rbound = rpre;
                        if(rpre == k){
                            rflag = false;
                        }
                    }
                }
            }
*/            
            for(int s = k + 1; s <= rbound; ++s){
                if(flag[s] == '0'){
                    flag[s] = '1';
                    maxFlag[s] = '1';
                    ++oneCount;
                }
            }
            // left
            if(rflag){
                lpVal = rpVal;
            }
            for(int t = k - 1; t >= 0; --t){
                if(winVec[t].pend >= xstart){
                    if(winVec[t].value >= lowLimit && winVec[t].value <= upLimit){
                        xstart = winVec[t].start - dx;
                        if(flag[t] == '0'){
                            lpre = lbound;
                            //maxFlag[t] = '1';
                            //++oneCount;
                            //
                            if(lpre < k){
                                float tval = wMean(winVec,t,lpre - 1,wmet);
                                float diff = tval - lpVal;
                                if(diff < 0){
                                    diff = -diff;
                                }
                                if(diff <= dy){
                                    lbound = t;
                                    //flag[t] = '1';
                                    //
                                    lflag = true;
                                }else{
                                    break;
                                }
                                
                            }else{
                                if(rflag){
                                    float tval = wMean(winVec,t,lpre - 1,wmet);
                                    float diff = tval - lpVal;
                                    if(diff < 0){
                                        diff = -diff;
                                    }
                                    if(diff <= dy){
                                        lbound = t;
                                        //flag[t] = '1';
                                        //
                                        lflag = true;
                                    }else{
                                        break;
                                    }
                                }else{
                                    if(winVec[lpre].start - winVec[t].pend > plateau){
                                        if(xlen >= plateau){
                                            // contain gap or another block
                                            break;
                                        }else{
                                            int ntLen = winVec[t].pend - winVec[t].start + 1;
                                            if(ntLen >= plateau){
                                                break;
                                            }else{
                                        
                                                if(lpre - t > 1){
                                                    int mid = (winVec[lpre].start + winVec[t].pend) / 2;
                                                    int rnum = 0,lnum = 0;
                                                    int rStart = 0;
                                                    for(int u = t + 1; u < lpre; ++u){
                                                        if(winVec[u].pend <= mid){
                                                            ++lnum;
                                                        }else if(winVec[u].start >= mid){
                                                            if(rnum == 0){
                                                                rStart = u;
                                                            }
                                                            ++rnum;
                                                        }
                                                    }
                                                    if(lnum == 0){
                                                        if(rnum == 0){
                                                            // cross middle
                                                            lbound = t;
                                                            //flag[t] = '1';
                                                            //
                                                            lpVal = wMean(winVec,t,lpre,wmet);
                                                            lflag = true;
                                                        }else{
                                                            break;
                                                        }
                                                    }else{
                                                        break;
                                                        /*
                                                        if(rnum == 0){
                                                            break;
                                                        }else{
                                                            float lw = wMean(winVec,t + 1,t + lnum);
                                                            float rw = wMean(winVec,rStart,lpre - 1);
                                                            float mdiff = rw - lw;
                                                            if(mdiff < 0){
                                                                mdiff = -mdiff;
                                                            }
                                                            if(mdiff <= dy){
                                                                lbound = t;
                                                                lpVal = wMean(winVec,t,lpre);
                                                                lflag = true;
                                                            }
                                                        }
                                                        */
                                                    }
                                                }else{
                                                    if(winVec[lpre].start - winVec[t].pend > plateau){
                                                        // gap
                                                        break;
                                                    }else{
                                                        lbound = t;
                                                        //flag[t] = '1';
                                                        //
                                                        lpVal = wMean(winVec,t,lpre,wmet);
                                                        lflag = true;
                                                    }
                                                }
                                            }
                                        }
                                    }else{
                                        float tval = wMean(winVec,t,lpre,wmet);
                                        if(xlen >= plateau){
                                            float diff = winVec[lpre].value - tval;
                                            if(diff < 0){
                                                diff = -diff;
                                            }
                                            if(diff <= dy){
                                                lbound = t;
                                                //flag[t] = '1';
                                                lpVal = tval;
                                                lflag = true;
                                            }else{
                                                break;
                                            }
                                        }else{
                                            int ntLen = winVec[t].pend - winVec[t].start + 1;
                                            if(ntLen >= plateau){
                                                float diff = winVec[t].value - tval;
                                                if(diff < 0){
                                                    diff = -diff;
                                                }
                                                if(diff <= dy){
                                                    lbound = t;
                                                    //flag[t] = '1';
                                                    lpVal = tval;
                                                    lflag = true;
                                                }else{
                                                    break;
                                                }
                                            }else{
                                                lbound = t;
                                                //flag[t] = '1';
                                                lpVal = tval;
                                                lflag = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }else{
                            //lmark = true;
                            break;
                        }
                    }else{
                        if(flag[t] == '1'){
                            //lmark = true;
                            break;
                        }
                    }
                }else{
                    //if(flag[t] == '1'){
                        //lmark = true;
                    //}
                    break;
                }
            }
            //
/*            
            if(lflag){
                if(lpre - lbound > 1){
                    if(winVec[lpre].start - winVec[lbound].pend > plateau){
                        if(winVec[lbound].pend - winVec[lbound].start < plateau){
                            flag[lbound] = '0';
                        }
                        lbound = lpre;
                        if(lpre == k){
                            lflag = false;
                        }
                    }
                }
            }
*/            
            for(int u = lbound; u < k; ++u){
                if(flag[u] == '0'){
                    flag[u] = '1';
                    maxFlag[u] = '1';
                    ++oneCount;
                }
            }
            //
            //char type = '0';
            if(xlen >= plateau || rflag || lflag){
                flag[k] = '1';
                MnodeRan tRan = {lbound,rbound};
                mran.push_back(tRan);
            }
        }else{
            break;
        }

    }
    //
    if(mran.empty()){
        lpos = 0;
        rpos = n - 1;
        //eName = "bottom";
        return;
    }
    sort(mran.begin(),mran.end(),posSort); 
//------------------------
    //
    int rSize = mran.size();
    
    bool outLast = false;
    if(mran.back().posEnd < n - 1){
        lpos = mran.back().posEnd + 1;
        rpos = n - 1;
        outLast = true;
    }else{
        lpos = mran.back().posStart;
        rpos = mran.back().posEnd;
        // not output
    }
    
    int prePos = 0;
    //string name = "";
    for(int r = 0; r < rSize; ++r){
        int lowStart = prePos;
        bool lowFlag = true;
        if(prePos < mran[r].posStart){
            for(int z = prePos; z <= mran[r].posStart - 1; ++z){
                if(winVec[z].value > baseValue){
                    if(lowFlag && z > lowStart){
                        storeValue(winVec,ass,tchr,lowStart,z - 1,wmet,ofh);
                    }
                    storeValue(winVec,ass,tchr,z,z,wmet,ofh);
                    lowFlag = false;
                }else{
                    if(! lowFlag){
                        lowStart = z;
                        lowFlag = true;
                    }
                }
            }
            if(lowFlag){
                storeValue(winVec,ass,tchr,lowStart,mran[r].posStart - 1,wmet,ofh);
            }
        }
        //
        
        if(r < rSize - 1){
            storeValue(winVec,ass,tchr,mran[r].posStart,mran[r].posEnd,wmet,ofh);
        }else{
            if(outLast){
                storeValue(winVec,ass,tchr,mran[r].posStart,mran[r].posEnd,wmet,ofh);
            }
        }
        //
        prePos = mran[r].posEnd + 1;
    }
}

void scanNode(string &inFile,int win,int xdis,float ydis,float baseValue,int plateau,int wmet,string &outFile){
    ifstream in(inFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<inFile<<endl;
        exit(1);
    }
    ofstream out(outFile.c_str());
    if(! out){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    string ass,tChr;
    int start,pend;
    float value;
    vector<WinDiv> winVec;
    int count = 0;
    string preChr = "",preAss = "";
    string vaLine;
    stringstream strStream;
    //string eName = "";
    int nnode;
    while(getline(in,vaLine)){
        strStream<<vaLine;
        strStream >> ass;
        strStream >> tChr;
        strStream >> start;
        strStream >> pend;
        strStream >> value;
        strStream >> nnode;
        
        strStream.clear();
        strStream.str("");
        //
        if(preAss != "" && ass != preAss){            
            int lpos = 0,rpos = 0;
            if(! winVec.empty()){
                mergeNode(winVec,preAss,preChr,xdis,ydis,baseValue,plateau,wmet,out,lpos,rpos);
                storeValue(winVec,preAss,preChr,lpos,rpos,wmet,out);
            }
            //
            winVec.clear();
            count = 0;
        }else{
            if(preChr != "" && tChr != preChr){
                int lpos = 0,rpos = 0;
                if(! winVec.empty()){
                    mergeNode(winVec,ass,preChr,xdis,ydis,baseValue,plateau,wmet,out,lpos,rpos);
                    storeValue(winVec,ass,preChr,lpos,rpos,wmet,out);
                }
                //
                winVec.clear();
                count = 0;
            }else{
                if(count == win){
                    int lpos = 0,rpos = 0;
                    mergeNode(winVec,ass,tChr,xdis,ydis,baseValue,plateau,wmet,out,lpos,rpos);
                    vector<WinDiv> tWin;
                    for(size_t x = lpos; x <= rpos; ++x){
                        tWin.push_back(winVec[x]);
                    }
                    //
                    winVec.clear();
                    for(size_t m = 0; m < tWin.size(); ++m){
                        winVec.push_back(tWin[m]);
                    }
                    count = 0;
                }
            }
        }
        //
        //if(value >= lowVal && value <= upVal){
            WinDiv tval = {start,pend,value,nnode};
            winVec.push_back(tval);
            ++count;
        //}
        preChr = tChr;
        //
        preAss = ass;
    }
    int lpos = 0,rpos = 0;
    if(! winVec.empty()){
        mergeNode(winVec,preAss,preChr,xdis,ydis,baseValue,plateau,wmet,out,lpos,rpos);
        storeValue(winVec,preAss,preChr,lpos,rpos,wmet,out);
    }
    in.close();
    out.close();
}

void mergeTop(char *inFile,int plateau,float topth,int wmet,string &outFile){
    ifstream in(inFile);
    if(! in){
        cerr<<"Error: file open failed. "<<inFile<<endl;
        exit(1);
    }
    ofstream out(outFile.c_str());
    if(! out){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    //
    string ass,tChr;
    int start,pend;
    float value;
    int nnode;
    
    vector<WinDiv> winVec;
    WinDiv firDiv;
    int count = 0;
    string preChr = "",preAss = "";
    string vaLine;
    stringstream strStream;

    //
    bool firFlag = false;
    
    int midSum = 0;
    while(getline(in,vaLine)){
        strStream<<vaLine;
        strStream >> ass;
        strStream >> tChr;
        strStream >> start;
        strStream >> pend;
        strStream >> value;
        strStream >> nnode;
        
        strStream.clear();
        strStream.str("");
        //
        int len = pend - start + 1;
        if(ass != preAss){
            if(firFlag){
                out<<preAss<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                for(auto &tdiv : winVec){
                    out<<preAss<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                }
                //
                winVec.clear();
            }
            //
            if(value < topth){
                out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                firFlag = false;
            }else{
                firDiv.start = start;
                firDiv.pend = pend;
                firDiv.value = value;
                firDiv.nnode = nnode;
                firFlag = true;
                midSum = 0;
            }
        }else{
            if(tChr != preChr){
                if(firFlag){
                    out<<ass<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                    for(auto &tdiv : winVec){
                        out<<ass<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                    }
                    //
                    winVec.clear();
                }
                //midSum = 0;
                if(value < topth){
                    out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                    firFlag = false;
                }else{
                    firDiv.start = start;
                    firDiv.pend = pend;
                    firDiv.value = value;
                    firDiv.nnode = nnode;
                    firFlag = true;
                    midSum = 0;
                }
                
            }else{
                if(! firFlag){
                    if(value < topth){
                        out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                    }else{
                        firDiv.start = start;
                        firDiv.pend = pend;
                        firDiv.value = value;
                        firDiv.nnode = nnode;
                        firFlag = true;
                        midSum = 0;
                    }
                }else{
                    if(value < topth){
                        midSum += len;
                        if(midSum < plateau){
                            WinDiv xdiv = {start,pend,value,nnode};
                            winVec.push_back(xdiv);
                        }else{
                            out<<ass<<"\t"<<tChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                            for(auto &tdiv : winVec){
                                out<<ass<<"\t"<<tChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                            }
                            out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                            //
                            firFlag = false;
                            winVec.clear();
                        }
                    }else{
                        if(start - firDiv.pend <= plateau){
                            int midLen = 0;
                            double mValue = 0;
                            int tlNode = 0;
                            float tv2;
                            if(wmet == 1){
                                int firLen = firDiv.pend - firDiv.start + 1;
                                for(auto &tdiv : winVec){
                                    int xlen = tdiv.pend - tdiv.start + 1;
                                    mValue += tdiv.value * xlen;
                                    midLen += xlen;
                                    tlNode += tdiv.nnode;
                                }
                                //
                                tv2 = ((double)firDiv.value * firLen + mValue + (double)value * len) / (firLen + midLen + len);
                            }else{
                                int firLen = firDiv.nnode;
                                for(auto &tdiv : winVec){
                                    int xlen = tdiv.nnode;
                                    mValue += tdiv.value * xlen;
                                    midLen += xlen;
                                    tlNode += tdiv.nnode;
                                }
                                //
                                tv2 = ((double)firDiv.value * firLen + mValue + (double)value * nnode) / (firLen + midLen + nnode);
                            }
                            firDiv.pend = pend;
                            firDiv.value = tv2;
                            firDiv.nnode += tlNode + nnode;
                            
                        }else{
                            out<<ass<<"\t"<<tChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                            for(auto &tdiv : winVec){
                                out<<ass<<"\t"<<tChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                            }
                            //
                            firDiv.start = start;
                            firDiv.pend = pend;
                            firDiv.value = value;
                            firDiv.nnode = nnode;
                        }
                        //
                        winVec.clear();
                        midSum = 0;
                    }
                }
            }
        }
        //
        preChr = tChr;
        preAss = ass;
    }
    //
    if(firFlag){
        out<<preAss<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
        for(auto &tdiv : winVec){
            out<<preAss<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
        }
        //
        winVec.clear();
    }
    //
    in.close();
    out.close();
}

// slop/plateau,peak,slop/plateau -> plateau
// first WinDiv, middle WinDid vector
// middle block <= 10
// middle block sum <= 10
// flank >= max(50,plateau / 5)
// flank sum / middle block >= 100

void finalNode(string &inFile,int xdis,float ydis,int plateau,int wmet,string &outFile){
    ifstream in(inFile.c_str());
    if(! in){
        cerr<<"Error: file open failed. "<<inFile<<endl;
        exit(1);
    }
    ofstream out(outFile.c_str());
    if(! out){
        cerr<<"Error: file open failed. "<<outFile<<endl;
        exit(1);
    }
    //
    string ass,tChr;
    int start,pend;
    float value;
    int nnode;
    
    vector<WinDiv> winVec;
    WinDiv firDiv;
    int count = 0;
    string preChr = "",preAss = "";
    string vaLine;
    stringstream strStream;

    //
    bool firFlag = false;
    
    int midSum = 0;
    while(getline(in,vaLine)){
        strStream<<vaLine;
        strStream >> ass;
        strStream >> tChr;
        strStream >> start;
        strStream >> pend;
        strStream >> value;
        strStream >> nnode;
        
        strStream.clear();
        strStream.str("");
        //
        int len = pend - start + 1;
        if(ass != preAss){
            if(firFlag){
                out<<preAss<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                for(auto &tdiv : winVec){
                    out<<preAss<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                }
                //
                winVec.clear();
            }
            //
            if(len < plateau){
                out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                firFlag = false;
            }else{
                firDiv.start = start;
                firDiv.pend = pend;
                firDiv.value = value;
                firDiv.nnode = nnode;
                firFlag = true;
                midSum = 0;
            }
        }else{
            if(tChr != preChr){
                if(firFlag){
                    out<<ass<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                    for(auto &tdiv : winVec){
                        out<<ass<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                    }
                    //
                    winVec.clear();
                }
                //midSum = 0;
                if(len < plateau){
                    out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                    firFlag = false;
                }else{
                    firDiv.start = start;
                    firDiv.pend = pend;
                    firDiv.value = value;
                    firDiv.nnode = nnode;
                    firFlag = true;
                    midSum = 0;
                }
                
            }else{
                if(! firFlag){
                    if(len < plateau){
                        out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                    }else{
                        firDiv.start = start;
                        firDiv.pend = pend;
                        firDiv.value = value;
                        firDiv.nnode = nnode;
                        firFlag = true;
                        midSum = 0;
                    }
                }else{
                    if(len < plateau){
                        midSum += len;
                        if(midSum < plateau){
                            WinDiv xdiv = {start,pend,value,nnode};
                            winVec.push_back(xdiv);
                        }else{
                            out<<ass<<"\t"<<tChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                            for(auto &tdiv : winVec){
                                out<<ass<<"\t"<<tChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                            }
                            out<<ass<<"\t"<<tChr<<"\t"<<start<<"\t"<<pend<<"\t"<<fixed<<setprecision(5)<<value<<"\t"<<nnode<<endl;
                            //
                            firFlag = false;
                            winVec.clear();
                        }
                    }else{
                        bool flag = false;
                        if(start - firDiv.pend <= xdis){
                            float diff = value - firDiv.value;
                            if(diff < 0){
                                diff = -diff;
                            }
                            if(diff <= ydis){
                                int midLen = 0;
                                double mValue = 0;
                                int tlNode = 0;
                                float tv1,tv2;
                                if(wmet == 1){
                                    int firLen = firDiv.pend - firDiv.start + 1;
                                    for(auto &tdiv : winVec){
                                        int xlen = tdiv.pend - tdiv.start + 1;
                                        mValue += tdiv.value * xlen;
                                        midLen += xlen;
                                        tlNode += tdiv.nnode;
                                    }
                                    //
                                    tv1 = ((double)firDiv.value * firLen + (double)value * len) / (firLen + len);
                                    tv2 = ((double)firDiv.value * firLen + mValue + (double)value * len) / (firLen + midLen + len);
                                }else{
                                    int firLen = firDiv.nnode;
                                    for(auto &tdiv : winVec){
                                        int xlen = tdiv.nnode;
                                        mValue += tdiv.value * xlen;
                                        midLen += xlen;
                                        tlNode += tdiv.nnode;
                                    }
                                    //
                                    tv1 = ((double)firDiv.value * firLen + (double)value * nnode) / (firLen + nnode);
                                    tv2 = ((double)firDiv.value * firLen + mValue + (double)value * nnode) / (firLen + midLen + nnode);
                                }
                                    
                                float df2 = tv1 - tv2;
                                if(df2 >= -ydis && df2 <= ydis){
                                    firDiv.pend = pend;
                                    firDiv.value = tv2;
                                    firDiv.nnode += tlNode + nnode;
                                }else{
                                    flag = true;
                                }
                            }else{
                                flag = true;
                            }
                        }else{
                            flag = true;
                        }
                        //
                        if(flag){
                            out<<ass<<"\t"<<tChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
                            for(auto &tdiv : winVec){
                                out<<ass<<"\t"<<tChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
                            }
                            //
                            firDiv.start = start;
                            firDiv.pend = pend;
                            firDiv.value = value;
                            firDiv.nnode = nnode;
                        }
                        //
                        winVec.clear();
                        midSum = 0;
                    }
                }
            }
        }
        //
        preChr = tChr;
        preAss = ass;
    }
    //
    if(firFlag){
        out<<preAss<<"\t"<<preChr<<"\t"<<firDiv.start<<"\t"<<firDiv.pend<<"\t"<<fixed<<setprecision(5)<<firDiv.value<<"\t"<<firDiv.nnode<<endl;
        for(auto &tdiv : winVec){
            out<<preAss<<"\t"<<preChr<<"\t"<<tdiv.start<<"\t"<<tdiv.pend<<"\t"<<fixed<<setprecision(5)<<tdiv.value<<"\t"<<tdiv.nnode<<endl;
        }
        //
        winVec.clear();
    }
    //
    in.close();
    out.close();
}

void getBase(string &inFile,float yper,float &baseVal,float &yDis,float &topth,float &maxAll){
    //
    ifstream in(inFile);
    string vaLine;
    stringstream strStream;
    string ass,tChr;
    int start,pend;
    float value;
    //
    vector<float> allValue;
    while(getline(in,vaLine)){
        strStream<<vaLine;
        strStream >> ass;
        strStream >> tChr;
        strStream >> start;
        strStream >> pend;
        strStream >> value;
        if(strStream.bad()){
            cerr<<"Error: please check the input file. "<<inFile<<endl;
            cerr<<"Error line: "<<vaLine<<endl;
            exit(1);
        }
        allValue.push_back(value);
        if(value > 1e7){
            cerr<<"Warning: value is too large, may be out of the resolution of the program. "<<vaLine<<endl;
        }
        //
        strStream.clear();
        strStream.str("");
    }
    //float minVal = 0.0f;
    in.close();
    //
    int upQuan = allValue.size() * 0.975 - 1;
    int downQuan = allValue.size() * 0.025 - 1;
    int upLoc = allValue.size() * 0.75 - 1;
    int downLoc = allValue.size() * 0.25 - 1;
    if(upQuan < 0){
        upQuan = 0;
    }
    if(downQuan < 0){
        downQuan = 0;
    }
    if(upLoc < 0){
        upLoc = 0;
    }
    if(downLoc < 0){
        downLoc = 0;
    }
    //sort(allValue.begin(),allValue.end());
    vector<float>::iterator it;
    float val25 = 0.0f,val75 = 0.0f,val975 = 0.0f; 
    
    nth_element(allValue.begin(),allValue.begin() + downQuan,allValue.end());
    it = allValue.begin() + downQuan;
    baseVal = *it;
/*    
    nth_element(allValue.begin(),allValue.begin() + downLoc,allValue.end());
    it = allValue.begin() + downLoc;
    val25 = *it;
    
    nth_element(allValue.begin(),allValue.begin() + upLoc,allValue.end());
    it = allValue.begin() + upLoc;
    val75 = *it;
*/
    nth_element(allValue.begin(),allValue.begin() + upQuan,allValue.end());
    it = allValue.begin() + upQuan;
    val975 = *it;
    
    yDis = (val975 - baseVal) * (yper / 100);
    //
    //float dis = val75 - val25;
    //upValue = val75 + 1.5 * dis;
    //downValue = val25 - 1.5 * dis;
    int topn = allValue.size() * 0.95 - 1;
    nth_element(allValue.begin(),allValue.begin() + topn,allValue.end());
    it = allValue.begin() + topn;
    topth = *it;
    maxAll = *max_element(allValue.begin(),allValue.end());
} 
  
// minimal distance between two adjacent blocks on the y-axis 
void merge_usage(){
    cout<<"Usage: BLOOM --input input_file --outdir output_dir"<<endl;
    cout<<"--input      File    input file"<<endl;
    cout<<"--outdir     Dir     output directory"<<endl;
    cout<<"--prefix     String  prefix of output file"<<endl;
    cout<<"--xdis       Int     maximal distance allowed when searching points with similar values, by default: 500"<<endl;
    cout<<"--ydis       Float   threshold value in decimal to seperate two adjacent blocks (0,), forbid --yper"<<endl;
    cout<<"--yper       Float   threshold value in percentage to seperate two adjacent blocks (0,100), by default: (0.975 quantile - 0.025 quantile) * 0.05, forbid --ydis"<<endl;
    cout<<"--bvalue     Float   threshold of bottom, by default: boundary value at lower 2.5%"<<endl;
    cout<<"--pwidth     Int     maximal distance on x-axis allowed when merging two adjacent points with similar values on y-axis directly, by default 50"<<endl;
    cout<<"--iteration  Int     times of iteration, by default 5"<<endl;
    //cout<<"--filter             filter out abnormal value before block seeking"<<endl;
    cout<<"--wmet       Int     methods to calculate weighted value. 1 -- by segment size, 2 -- by number of data points"<<endl; 
    cout<<"--mtop               merge top 5%"<<endl;
    cout<<"--help               "<<endl;
}
//inFile,xdis,ydis,baseValue,outFile
int main(int argc,char **argv){
    char *inFile = nullptr,*outDir = nullptr;
    string prefix = "merge";
    float ydis = 0.05f, bvalue = 0.05f;
    float yper = 5.0f;
    int xdis = 500;
    int pwidth = 50, iter = 5;
    int wmet = 1;
    
    bool bFlag = false,yFlag = false;
    bool mtop = false;
    
    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i],"--input") == 0 || strcmp(argv[i],"--input") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --input is used but there is no value"<<endl;
                return 1;
            }
            inFile = argv[i];
        }else if(strcmp(argv[i],"--outdir") == 0 || strcmp(argv[i],"--outdir") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --outdir is used but there is no value"<<endl;
                return 1;
            }
            outDir = argv[i];
        }else if(strcmp(argv[i],"--prefix") == 0 || strcmp(argv[i],"--prefix") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --prefix is used but there is no value"<<endl;
                return 1;
            }
            prefix = argv[i];
        }else if(strcmp(argv[i],"--xdis") == 0 || strcmp(argv[i],"--xdis") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --xdis is used but there is no value"<<endl;
                return 1;
            }
            xdis = atoi(argv[i]);
        }else if(strcmp(argv[i],"--ydis") == 0 || strcmp(argv[i],"--ydis") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --ydis is used but there is no value"<<endl;
                return 1;
            }
            ydis = atof(argv[i]);
            yFlag = true;
        }else if(strcmp(argv[i],"--yper") == 0 || strcmp(argv[i],"--yper") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --yper is used but there is no value"<<endl;
                return 1;
            }
            yper = atof(argv[i]);
        }else if(strcmp(argv[i],"--bvalue") == 0 || strcmp(argv[i],"--bvalue") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --bvalue is used but there is no value"<<endl;
                return 1;
            }
            bvalue = atof(argv[i]);
            bFlag = true;
        }else if(strcmp(argv[i],"--pwidth") == 0 || strcmp(argv[i],"--pwidth") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --pwidth is used but there is no value"<<endl;
                return 1;
            }
            pwidth = atoi(argv[i]);
        }else if(strcmp(argv[i],"--iteration") == 0 || strcmp(argv[i],"--iteration") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --iteration is used but there is no value"<<endl;
                return 1;
            }
            iter = atoi(argv[i]);
        }else if(strcmp(argv[i],"--wmet") == 0 || strcmp(argv[i],"--wmet") == 0){
            ++i;
            if(i == argc){
                cerr<<"Error: --wmet is used but there is no value"<<endl;
                return 1;
            }
            wmet = atoi(argv[i]);
        }else if(strcmp(argv[i],"--mtop") == 0 || strcmp(argv[i],"-mtop") == 0){
            mtop = true;
        }else if(strcmp(argv[i],"--mtop") == 0 || strcmp(argv[i],"-mtop") == 0){
            merge_usage();
            return 1;
        }else{
            cerr<<"Error: undefined parameter"<<endl;
            merge_usage();
            return 1;
        }
    }
    //
    if(inFile == nullptr){
        cerr<<"Error: --input is required"<<endl;
        exit(1);
    }
    if(outDir == nullptr){
        cerr<<"Error: --outDir is required"<<endl;
        exit(1);
    }
    if(! opendir(outDir)){
        if(mkdir(outDir,S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)){
            cerr<<"Error: failed to create output directory"<<endl;
            exit(1);
        }
    }
    string preValFile = inFile;
    string tdir = outDir;
    //
    float botVal = 0.0f, ryDis = 0.0f;
    //upVal = 0.0f, downVal = 0.0f
    float topth = 0.0f,maxAll = 0.0f;
    getBase(preValFile,yper,botVal,ryDis,topth,maxAll);
    if(! bFlag){
        bvalue = botVal;
    }
    
    if(! yFlag){
        ydis = ryDis;
    }
    //
    //if(ydis < 1e-5){
    //    cerr<<"Warning: value was too small. --ydis was set to 1e-5"<<endl;
    //    ydis = 1e-5;
    //}
    //cout<<"upper quartile : "<<upVal<<endl;
    //cout<<"lower quartile : "<<downVal<<endl;
    cout<<"Max value : "<<maxAll<<endl;
    cout<<"Threshold at top 5% : "<<topth<<endl;
    cout<<"--bvalue : "<<bvalue<<endl;
    if(! bFlag){
        cout<<"--yper : "<<yper<<endl;
    }
    cout<<"--ydis : "<<ydis<<endl;
    cout<<"--xdis : "<<xdis<<endl;
    cout<<"--pwidth : "<<pwidth<<endl;
    cout<<"--wmet : "<<wmet<<endl;
    cout<<"------------------"<<endl;
    
    if(mtop){
        preValFile = tdir + "/" + prefix + ".mtop.out";
        mergeTop(inFile,pwidth,topth,wmet,preValFile);
    }
    //
    int preCount = 0;
    int k = 0;
    for(k = 0; k < iter;++k){
        cout<<"iteration "<<k+1<<endl;
        saveCount = 0;
        string outFile = tdir + "/" + prefix + ".iter." + to_string(k + 1) + ".out";
        scanNode(preValFile,100,xdis,ydis,bvalue,pwidth,wmet,outFile);
        
        preValFile = outFile;
        //
        if(saveCount == preCount){
            break;
        }
        preCount = saveCount;
    }
    string fnFile = tdir + "/" + prefix + ".final." + to_string(k+1) + ".out";
    finalNode(preValFile,xdis,ydis,pwidth,wmet,fnFile);
}






    


