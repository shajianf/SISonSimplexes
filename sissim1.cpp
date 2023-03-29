//first test for sis on simplex. 7/11/1401 (run is good! but not sure for result)
//add true way for generate simplex. 29.11.1401
//finally successful! 6.12.1401
#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<cstdlib>
#include<vector>

#define N 10000
#define M 5

using namespace std;

ofstream Time;
ofstream Pk;
vector <vector<int>> Sim22(N);              //Adjecancy matrix for 2-simplexes..
vector <vector<int>> Sim1_up(N);            //Adjecancy matrix for 1-simplexes.

double k = 20;                              //degree of 1-simplexes.
double k3 = 6;                              //degree of 2-simplexes.

double p1 = ((k - (2*k3))/((N-1)-(2*k3)));              //probability of create 1-simplex.
double p2 = (2*k3)/((N-1)*(N-2));                       //probability of create 2-simplex.


int N0 = 100;                       // Initial infected nodes..

int ensemble = 10;
int tmax = 1000;                     // time steps.

double beta = 0.;                    // probability of a node get infected (for !-simplex)
double beta3 = 0.026;                 //probability of a node get infected (for 2-simplexes).
double muo = 0.20;                    //probability of recovery.

int sis[N]={0};                     //state of each node; 1= infected , 0= suseptible.at time t. 
int sis_up[N]={0};                  //state of each node at time t+1.

int inf_sim1[N]={0};                //number of infected neighbors (for 1-simplexes).
int inf_sim2[N]={0};                //number of infected neighbors (for 2-simplexes).

void add(int a0, int b0,int c0){    // for add 1 link between b & c.

    bool flg0 = false;
    for (int s=0; s< Sim1_up[b0].size();s++){

        int q01=0;
        q01 = Sim1_up[b0][s];

        
        if ( q01 == c0) {
            flg0 = false;
            break;
        }
    }

    if (flg0){
            
        Sim1_up[c0].push_back(b0);              //add to 1-simplexes.
        Sim1_up[b0].push_back(c0);


        Sim22[a0].push_back(b0);                //add to 2-simplexes.
        Sim22[a0].push_back(c0);

        Sim22[b0].push_back(a0);
        Sim22[b0].push_back(c0);

        Sim22[c0].push_back(a0);
        Sim22[c0].push_back(b0);
        
        //Pk << x <<"\t" << a0 << "\t" << b0 << "\t" << c0 <<endl;
    }

    
}

void add2 (int a1,int b1, int c1){              //for add 2 links between a & c and b& c.


    bool flg = true;
    bool flg1 = true;
    for (int s=0; s< Sim1_up[a1].size();s++){

        int q01=0;
        q01 = Sim1_up[a1][s];

        if ( q01 == c1) {
            flg = false;
            break;
        }


    }

    for (int s1=0; s1< Sim1_up[b1].size();s1++){

        int q01=0;
        q01 = Sim1_up[b1][s1];

        
        if ( q01 == c1) {
            flg1 = false;
            break;
        }
    }

    if ((flg)&&(flg1)){

        Sim1_up[a1].push_back(c1);
        Sim1_up[c1].push_back(a1);

        Sim1_up[b1].push_back(c1);
        Sim1_up[c1].push_back(b1);



        Sim22[a1].push_back(b1);
        Sim22[a1].push_back(c1);

        Sim22[b1].push_back(a1);
        Sim22[b1].push_back(c1);

        Sim22[c1].push_back(a1);
        Sim22[c1].push_back(b1);
 
    }

}

void add3 (int a2,int b2, int c2){              //add 3 links between a & b & c.


    bool flg = true;
    bool flg1 = true;
    bool flg2 = true;
    for (int s=0; s< Sim1_up[a2].size();s++){

        int q01=0;
        q01 = Sim1_up[a2][s];

        if ( q01 == c2) {
            flg = false;
            break;
        }


    }

    for (int s1=0; s1< Sim1_up[b2].size();s1++){

        int q01=0;
        q01 = Sim1_up[b2][s1];

        
        if ( q01 == c2) {
            flg1 = false;
            break;
        }
    }

    for (int s2=0; s2< Sim1_up[b2].size();s2++){

        int q01=0;
        q01 = Sim1_up[b2][s2];

        
        if ( q01 == a2) {
            flg2 = false;
            break;
        }
    }

    if ((flg)&&(flg1)&&(flg2)){

        Sim1_up[a2].push_back(c2);
        Sim1_up[c2].push_back(a2);

        Sim1_up[b2].push_back(c2);
        Sim1_up[c2].push_back(b2);

        Sim1_up[a2].push_back(b2);
        Sim1_up[b2].push_back(a2);


        Sim22[a2].push_back(b2);
        Sim22[a2].push_back(c2);

        Sim22[b2].push_back(a2);
        Sim22[b2].push_back(c2);

        Sim22[c2].push_back(a2);
        Sim22[c2].push_back(b2);

    }

}


int main(){

    clock_t runtime = clock();
    Time.open("SISonSimtestS0.8.txt");
    //Pk.open("testbeta.txt");
    srand(time(NULL));

    while ( beta < 0.026 ){
        //cout<<"++++++++++++++++++"<< " " << beta <<endl;

        int ens[ensemble] ={0};
     
        for(int e=0 ; e <ensemble ; e++){  

            for(int b =0; b < N; b++){          //reset

                Sim1_up[b].clear();
                Sim22[b].clear();
                sis [b]=0;
                sis_up [b]=0;

            }

            Sim1_up.resize(N);
            Sim22.resize(N);
            //generate Simplexes
            //-------------------------------- //generate 1-simplexes by p1 as an ER network.
            for(int i=0; i<N ; i++){            
                double num_edge =0;
                for(int j=0; j<i; j++){


                    float ran1 = rand()/(1.+RAND_MAX);
            
                    if(ran1 < p1 ){

                        Sim1_up[j].push_back(i);
                        Sim1_up[i].push_back(j);

                    }
                }
       
            }

            // ------------------------ //generate 1-simplexes by <k3>.
            double k3_test =0;

            while (k3_test < k3) {         

                int i = 0;              //choosing 3 nodes randomly.
                int j = 0; 
                int k = 0;

                i = rand()%N;
                j = rand()%N;
                k = rand()%N;

                if ((i != j)&&(j != k) && (k != i)){            // Not SELFEDGES.

                    bool flag1 = false;         // i,j.

                    for (int s0 =0; s0< Sim1_up[i].size(); s0++){       //find link between i & j.
                        int q00 = 0;
                        q00 = Sim1_up[i][s0];

                        if (q00 == j){
                    
                            flag1 = true;
                            break;
                        }
                    }
        
                    bool flag2 = false;         // i,k.

                    for(int s=0; s< Sim1_up[i].size();s++){         //find link between i & k.

                        int q10=0;
                        q10 = Sim1_up[i][s];
                        
                        if (q10==k){

                    
                            flag2 = true;
                            break;
                        }
                        
                    }
            
                    bool flag3 = false;         //j,k

            
                    for (int s1=0; s1< Sim1_up[j].size(); s1++){    //find link between j & k.
                        int q20=0;

                        q20 = Sim1_up[j][s1];

                        if (q20==k){

                            flag3 = true;
                            break;
                        }

                    }

                    if ((flag1 == false) && (flag2 == false) && (flag3 == false)){      //there are not any link between i,j,k.
                
                        add3(i,j,k);            


                    }else if ((flag1)&& (flag2) && (flag3 == false)){       //there are 2 links between i,j & i,k.
                
                        add(i,j,k);
                        
                    }else if ((flag1)&& (flag2 == false) && (flag3)){
                
                        add(j,i,k);

                    }else if ((flag1 == false)&& (flag2) && (flag3)){
                
                        add(k,i,j);

                    }else if ((flag1 == false) && (flag2 == false) && (flag3)){ //there are 1 link between j,k.
                
                        add2(k,j,i);

                    }else if ((flag1 == false) && (flag2) && (flag3 == false)){
                
                        add2(i,k,j);

                    }else if ((flag1)&& (flag2==false) && (flag3==false)){
                    
                        add2(i,j,k);  
                    
                    }else{              //there are 3 links between i,j,k.

                        Sim22[i].push_back(j);
                        Sim22[i].push_back(k);

                        Sim22[j].push_back(i);
                        Sim22[j].push_back(k);

                        Sim22[k].push_back(i);
                        Sim22[k].push_back(j);
                    }

                    double kk=0;
                    for (int w =0; w < Sim22.size(); w++){      //cheking <k3> at any stps.

                        kk += (Sim22[w].size())/2.0;    
                    }
                    k3_test = kk/(N);
                    //cout << k3_test <<endl;
                }
            }//for while k3.
            
            // ------------------------------ end simplexes.

            cout << "Be happy!" <<endl;

            int y=0;
            while(y<N0 ){                       //initialization
                int O=0;
                O = rand()%N;
                if (sis_up[O]==0){
                    sis[O]=1;
                    //sis_up[O]=1;
                    y++;
                }
            } 

            int t=0;  
            while(t<tmax){

                for (int f=0; f<N; f++){
                    float ran1= rand()/(1.+RAND_MAX);

                    if ((sis[f]) && (ran1 < muo)){              //recovery with mu prob.
                        
                        sis_up[f]=0;
                    }  
                    
                    for(int h=0;h< Sim1_up[f].size(); h++){

                        bool flag =false;
                        int q=0;
                        q = Sim1_up[f][h];

                        float ran2 = rand()/(1.+RAND_MAX);

                        if((sis[f]==0) && (sis[q]==1)&&(ran2 <beta)){    //counting infected 1-simplexes's neighbors).

                            //inf_sim1[f]++;
                            sis_up[f] =1;
                            break;
                        }  
                     
                    }
                    if (sis_up[f]==0){
                        for (int c=0; c< Sim22[f].size();c= c+2){   //counting infected 2-simplexes's neighbors).
                    
                            int d = 0;
                            d = Sim22[f][c];

                            int d1 = 0;
                            d1 = Sim22[f][c+1];
                            float ran3 = rand()/(1.+RAND_MAX);

                            if ((sis[d]==1)&&(sis[d1]==1)&&(sis[f]==0)&&(ran3 <beta3)){

                                //inf_sim2[f] ++;
                                sis_up[f]=1;
                                break;

                            }
                        
            
                        }
                    }
                    /*int v =0;
                    int v1 =0;
                    
                    v = inf_sim1[f];        //1-simplexes infected neighbors.
                    v1 = inf_sim2[f];       //2-simplexes infected neighbors.

                    double cul = 1-((pow((1-beta),v))*(pow((1-beta3),v1))); // 1-(1-beta)^v*(1-beta3)^v1.
                    double ran2 = rand()/(1.+RAND_MAX);

                    if (ran2 < cul){
                        
                        sis_up[f] =1;
                        //Pk << v <<'\t' << v1 << "\t" << Sim1_up[f].size() <<endl;
                    }*/

                }//end f.

                //cout << "end" << endl;
                double var =0;
                for( int u=0; u<N ; u++){                               //change step.
                    
                    sis[u] = sis_up[u];
                    inf_sim1[u]=0;
                    inf_sim2[u]=0;
   
                }

                t++;

            }//end of while t.

            int densityofinf =0;
            for(int b=0; b < N; b++){

                densityofinf += sis_up[b];
            }

            ens[e] = densityofinf;
            
        }//end of ensemble.

        double inf =0;
        for(int a=0; a<ensemble ;a++){
            inf +=ens[a];
            ens[a]=0;
        }

        Time <<(beta*k)/muo <<"\t" << (beta) << '\t' << inf/double(N * ensemble )<< endl;
        cout << (beta*k)/muo <<"\t" << (beta) << '\t' << inf/double(N * ensemble )<< endl;

        beta += 0.001;
    }

    cout << "run time = " << (double) (clock() - runtime) / (CLOCKS_PER_SEC * 60) << " min" << endl;

 return 0;

}