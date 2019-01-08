#include <iostream>
#include <cstdio>
#include <time.h> 
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;



class ErlangDistribution {
	double shape;
	double lamda;
	public:
    ErlangDistribution ();
	ErlangDistribution (double,double);
    int getShape(){
    	return shape;
	}
	int getLamda(){
    	return lamda;
	}
	
};

class Polynomial{
	double coefficeint;
	int t_degree;
	public:
	Polynomial(double,int);
	double getCoefficient(){
		return coefficeint;
	}
	int getDegree(){
		return t_degree;
	}
	void setCoefficient(double a){
		coefficeint = a;
	}
	void setDegree(int b){
		t_degree = b;
	}
};

ErlangDistribution::ErlangDistribution () {
	//time_t t;
	//srand((unsigned) time(&t)+rand());
}

ErlangDistribution::ErlangDistribution (double a, double b) {
	shape = a; // ni
	lamda = b;
	time_t t;
	srand((unsigned) time(&t)+rand());
}

Polynomial::Polynomial(double a,int b){
	coefficeint = a;
	t_degree = b;
}

class Random{
    private:
    double generateRandomNumber(){
        return ((double) rand() / (RAND_MAX));
    }
    double calculatePI(int shape){
			double U = 1;
        	for(int i=0;i<shape;i++){
                double x = generateRandomNumber();
                if(x!=0){
                    U*=x;
                }
                else{
                    do{
                        x = generateRandomNumber();;
                    }while(x==0);
                    U*=x;
                }
            	
        	}
           // cout << "shape="<<shape<<" U="<< U<<endl;
        	return U;
	}
    public:
    Random();
    int chooseAPRandom(double arr[],int numAp){
        int index = -1;
        double rm = generateRandomNumber();
       // cout << "rm=" <<rm<<endl;
        if(rm<arr[0])
        {
            index = 0;
        }
        else{
            for(int i=1;i<numAp;i++)
            {
                if(rm>=arr[i-1] && rm<=arr[i])
                {
                    index = i;
                }
           
            }
        }
        return index;
    }
    int poissonRandomNumber(double lamda){
        int X=-1;
        double U =generateRandomNumber();
        double p=exp(-1*lamda), f=p;
        int i=0;
        do{
    
            if(U<f){
                //x=i;
                break;
            }
            p = double((lamda*p)/(i+1));
            f = f+p;
            i++;
        }while(U>=f);
        X=i;
        return X;
    }
    double erlangRandomnumber(double lamda,int n){
       // cout<< "n="<<n<<endl;
        double U = calculatePI(n);
        double rm = (-1/lamda)*log(U);
       // cout << "U="<<U <<" rm="<<rm <<endl;
        return rm;
    }
};
Random::Random(){
    time_t t;
	srand((unsigned) time(&t)+rand());
}

int FindMinIndex(double list[],int size){
    double min = list[0];
    int index = 0;
    for(int i=0;i<size;i++)
    {
        if(min>list[i]){
            min = list[i];
            index = i;
        }
    }
    return index;
}

double fractorial(int n){

	double  frac = 1;
    if(n==0||n==1)
    {
        frac = 1;
    }
    else
    {
        for(int i=1;i<=n;i++)
        {
            frac*=i;
        }
    }
	
	return frac;
}


double getProbValue(ErlangDistribution lists[],int i,int k,int size){

	//cout << "welcome to new model" << endl;
	double pi=1;

    if(size>0){
        
        // detemine which APs will be considered
        int indexCj[size-1]  ={0};
		int index = 0;
		for(int x=0;x<size;x++){
			if(i!=x)
			{
				//cout << "helloooo" <<endl;
				indexCj[index] = x;
				index++;
			}			
		}

        //cout << "find P[AP"<<i<<"]=min("<<"t"<<indexCj[0]<<",t"<<indexCj[1]<<")"<<endl;

        // initial all APs 
        ErlangDistribution er[size];

        er[0] = lists[i];
		double sumlamda =er[0].getLamda();
        
        for(int n=0;n<index;n++)
        {
            er[n+1] = lists[indexCj[n]];
			sumlamda += er[n+1].getLamda();
        }

        /*for(int x=0;x<size;x++)
        {
            cout << "n="<<er[x].getShape()<< " lamda=" << er[x].getLamda() <<endl;
        }*/

        //generate Cj ......
        int upperbound = 0;
        for(int i=1;i<size;i++){
           // cout << "k =" << k << " n="<< er[i].getShape() << endl;
           // upperbound *= ((k*er[i].getShape())-1);
		   upperbound+=((k*er[i].getShape())-1);
        }
		upperbound = upperbound+1;
       // cout << "there are Cj="<<upperbound <<endl;

		vector<vector<Polynomial>> poly;

		// collect each h(t)
		for(int n=1;n<size;n++){
			int nKn = ((k*er[n].getShape())-1);
			vector<Polynomial> arr;
			for(int j=0;j<=nKn;j++){
				double res = (double)pow(er[n].getLamda(),j)/(double)fractorial(j);
				Polynomial pl(res,j);
				arr.push_back(pl);
			}
			poly.push_back(arr);
		}

		//print
		/*cout << "poly size = "<< poly.size() << endl;
		for(int i=0;i<poly.size();i++){
			cout << "the " << (i+1) << "th "<<endl;
			cout << "=============================" <<endl;
			vector<Polynomial> arr = poly[i];
			cout << "arr size ="<< arr.size() << endl;
			for(int j=0;j<arr.size();j++)
			{
				Polynomial pl = arr[j];
				cout << "arr["<<j<<"]="<<pl.getCoefficient() << ", degree: t^" << pl.getDegree() <<endl;
			}
			cout << "=============================" <<endl;
		}*/

		// initial some cj
		vector<Polynomial> cj; 
	//	cout << "ddd" << endl;
		int indicateCj[upperbound];
		fill_n(indicateCj,upperbound,-1);
	//	cout << "hhhhh" << endl;
		vector<Polynomial> firstpoly = poly[0];
		for(int i=0;i<firstpoly.size();i++)
		{
			cj.push_back(firstpoly[i]);
			indicateCj[firstpoly[i].getDegree()] = i;
		}
		//cout << "first cj.size =" << cj.size() << ",,"<<upperbound <<endl;

		for(int i=1;i<poly.size();i++)
		{	
			//cout << "i="<<i <<endl;
			vector<Polynomial> Mj; 
			vector<Polynomial> nextpoly = poly[i];
			int indicateMj[upperbound];
			fill_n(indicateMj,upperbound,-1);
			
			int ti=1;
			int indexMj = 0;
			for(int x=0;x<cj.size();x++)
			{
		
				double ci = 0;
				
				for(int j=0;j<nextpoly.size();j++)
				{
	
					ci = (cj[x].getCoefficient()*nextpoly[j].getCoefficient());
					ti = cj[x].getDegree()+nextpoly[j].getDegree();
					//cout << "ti="<<ti<<endl;
					Polynomial poi(ci,ti);

					// store in mj
					if(indicateMj[ti]==-1) // no index in desried t_degree
					{
						indicateMj[ti]=indexMj;
						Mj.push_back(poi);
						indexMj++;
					}
					else{
						Mj[indicateMj[ti]].setCoefficient(Mj[indicateMj[ti]].getCoefficient()+poi.getCoefficient());
					}
				}
			}

			/*for(int k=0;k<Mj.size();k++)
			{
				cout << "Mj["<<k<<"]="<<Mj[k].getCoefficient() << " degree:"<< Mj[k].getDegree() << endl ;
			}*/	

			//update cj by mj
			int insertCj = cj.size();
			//cout << "Mj.size =" << Mj.size() << endl;
			for(int m=0;m<Mj.size();m++)
			{
			
				if(indicateCj[Mj[m].getDegree()]==-1)
				{
					cj.push_back(Mj[m]);
					indicateCj[Mj[m].getDegree()]=insertCj;
					insertCj++;
				}
				else
				{
					//update

					cj[indicateCj[Mj[m].getDegree()]].setCoefficient(Mj[m].getCoefficient());

					//cj[indicateCj[Mj[m].getDegree()]].setCoefficient(Mj[indicateCj[Mj[m].getDegree()]].getCoefficient());
				//	cout<< "hhhh:" << cj[indicateCj[Mj[m].getDegree()]].getCoefficient()<<endl;
				}
				//cout<< "kkkk"<<endl;
			}

		}

		//cout << "last cj.size =" << cj.size() <<endl;
		
		
		double p=0;
		for(int j=0;j<cj.size();j++)
		{
			double c21 = fractorial((k*er[0].getShape())+j-1);
			double c22 = pow(sumlamda,(k*er[0].getShape())+j);
		//	double c2=((double)(c[j]/fractorial(j)))*c21;
			double c2=double((cj[j].getCoefficient()*c21)/c22);
			//cout << "sumlamda="<< sumlamda <<" ,c2=" << c2 <<endl;
			//cout <<"c"<<j<<"="<<cj[j].getCoefficient()<< " (k"<<i<<"+j-1)!=("<<(k*er[0].getShape())+j-1<<")!="<<c21<< " cj*(k+j-1)!="<<c2<<endl;
			p+=c2;
		}
		//cout << "h(t)="<<p<<endl;
		double o11 = pow(er[0].getLamda(),(k*er[0].getShape()));
		double o12 = fractorial((k*(er[0].getShape()))-1);
		double o1 = o11/o12;
		//cout << "lamda"<<i<<"^k"<<i<<"="<<eri.getLamda()<<"^"<<eri.getShape()<<"="<<o11<<" (k"<<i<<"-1)!="<<"("<<((eri.getShape())-1)<<")!="<<o12<<" (lamda"<<i<<"^k"<<i<<")/(k"<<i<<"-1)!="<<o1<<endl;
		//cout << "o1=" << o1 << endl;
		pi*=(o1*p);

    }
    else{
        pi=-1;
    }


	return pi;
}

int main(int argc, char *argv[]) {
    int num_ap = 4;
    int N_Simulation = 1000;
    double file_size = 100; // Mbyte
    double bandwidth[num_ap] = {10,50,250,500};
    double count1[num_ap] = {0};
    double count2[num_ap] = {0};
    double ti1[num_ap]={0};
    double ti2[num_ap]={0};
    //cout << "hello new project" << endl;
    cout << "file_size = " << file_size << "Mbyte" <<endl;
    cout << "Bandwidth : " ;
    for(int i=0;i<num_ap;i++)
    {
        cout << "B["<<i<<"]=" <<bandwidth[i] << "MB/s ";
    }
    cout << endl;
    Random rm;

    double popsum  = 1;

    double arr[num_ap]={0};
    double lamda[num_ap]={0};
    double ni[num_ap]={0};
    double Eni[num_ap]={0};
    // generate poisson random number...
    
    cout << "************* Settings *************" << endl;
    cout << "Simulation Rounds :" << N_Simulation << endl;
    for(int n=0;n<N_Simulation;n++)
    {
    cout <<" sim round: "<<(n+1) <<"-th" << endl;
    //cout << "# of people used : ";
    for(int i=0;i<num_ap;i++)
    {
        lamda[i] = (double)(bandwidth[i]/file_size);
        //lamda[i] = (double)(file_size/bandwidth[i]);
        ni[i] = rm.poissonRandomNumber(lamda[i]);
        Eni[i]+=ni[i];
       // cout <<" AP"<<i<<"="<<ni[i];
    }
    //cout << endl;
        //rm.getListProbability(arr,num_ap);

        // without wifi selection
       // cout << "************** Without WIFI Selection **************" <<endl; 
        arr[0]= double(popsum/num_ap);
       // cout << "arr[0]="<<arr[0]<<endl;
        for(int i=1;i<num_ap;i++)
        {
            arr[i]=arr[i-1]+ double(popsum/num_ap);
         //cout << "arr["<<i<<"]="<<arr[i]<<endl;
        }

        int index  = rm.chooseAPRandom(arr,num_ap);
       // cout << "Select WIFI AP = " << index << endl;
        count1[index]+=1;
        //generate ti
        //cout << lamda[index] <<endl;
        double ti = rm.erlangRandomnumber(lamda[index],ni[index]+1);
        //cout << "download time ti = " << ti << endl;
        ti1[index]+=ti;
       //cout << "************** With WIFI Selection **************" <<endl; 
        //calculate min E[Ti]
        //cout << "Expected download time =";
        double eti[num_ap]={0};
        for(int i=0;i<num_ap;i++)
        {
            eti[i]=(ni[i]+1)/lamda[i];
            //cout << " eti["<< i << "]="<<eti[i];
        }
       // cout << endl;
        int selectwifi = FindMinIndex(eti,num_ap);
       // cout << "Select WIFI AP = " << selectwifi << endl;
        count2[selectwifi]+=1;
        double sti = rm.erlangRandomnumber(lamda[selectwifi],ni[selectwifi]+1);
       // cout << "download time ti = " << sti << endl;
        ti2[selectwifi]+=sti;
    }
    cout << "***********Simulation result****************"<<endl;
    cout << "----- Average # of user E[ni] -----" <<endl;
    double eeNi[num_ap]={0};
    for(int i=0;i<num_ap;i++)
    {
        double Ani = (double)(Eni[i]/N_Simulation);
        eeNi[i]=Ani;
        cout << "E[n"<<i<<"]="<<Ani << " ";
    }
    cout <<endl;
    cout << "==== without selection ====" <<endl;
    cout << "# of selected WiFi AP ";
    for(int i=0;i<num_ap;i++)
    {
        cout <<"AP["<<i<<"]="<<count1[i]<<endl;
    }

    cout << "Propability : " << endl;
    for(int i=0;i<num_ap;i++)
    {
        double p = (double)(count1[i]/N_Simulation);
        double eti = (double)(ti1[i]/N_Simulation);
        cout <<"P[AP"<<i<<"]="<<p<< " E[t"<<i<<"]="<<eti<<endl;
    }

    cout << "==== with selection ====" <<endl;
    cout << "# of selected WiFi AP ";
    for(int i=0;i<num_ap;i++)
    {
        cout <<"Select AP["<<i<<"]="<<count2[i]<<endl;
    }
    for(int i=0;i<num_ap;i++)
    {
        double p = (double)(count2[i]/N_Simulation);
        double eti = (double)(ti2[i]/N_Simulation);
        cout <<"P[AP"<<i<<"]="<<p<< " E[t"<<i<<"]="<<eti<<endl;
    }

    cout << "***********Math result****************"<<endl;

    ErlangDistribution es[num_ap]; // generate random variable

	for(int i=0;i<num_ap;i++){
		ErlangDistribution e (eeNi[i],lamda[i]);
		es[i]= e;
	}

    double Mprob[num_ap]={0};
	double p2 = 0,psum =0;
    for(int i=0;i<num_ap;i++){
		// new function
		p2 = getProbValue(es, i,1,num_ap);

        Mprob[i]=p2;
        cout << "P[AP"<<i<<" is the best]="<<p2<<"\n";
        psum+=p2;

    }

    return 0;
}