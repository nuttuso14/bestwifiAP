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
	int shape;
	double lamda;
	public:
    ErlangDistribution ();
	ErlangDistribution (int,double);
    int getShape(){
    	return shape;
	}
	double getLamda(){
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

ErlangDistribution::ErlangDistribution (int a, double b) {
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
				indexCj[index] = x;
				index++;
			}			
		}

        // initial all APs 
        ErlangDistribution er[size];

        er[0] = lists[i];
		double sumlamda =er[0].getLamda();
        
        for(int n=0;n<index;n++)
        {
            er[n+1] = lists[indexCj[n]];
			sumlamda += er[n+1].getLamda();
        }

        //generate Cj ......
        int upperbound = 0;
        for(int i=1;i<size;i++){
           // upperbound *= ((k*er[i].getShape())-1);
		   upperbound+=((k*er[i].getShape())-1);
        }
		upperbound = upperbound+1;

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
			double c2=double((cj[j].getCoefficient()*c21)/c22);
			p+=c2;
		}
		//cout << "h(t)="<<p<<endl;
		double o11 = pow(er[0].getLamda(),(k*er[0].getShape()));
		double o12 = fractorial((k*(er[0].getShape()))-1);
		double o1 = o11/o12;
		pi*=(o1*p);

    }
    else{
        pi=-1;
    }


	return pi;
}


int main(int argc, char *argv[]) {

    if(argc<8){
        cerr << "Usage: " << argv[0] << " <SIM_ROUND> " <<" <Num_AP>" <<" <File_size>"  << " <Bandwidth1> " << " <Bandwidth2> " << " ... "<< " <Lambda1> " << " <Lamba2> " << " ... " << endl;
		return 1;
    }

   
    int N_Simulation = atoi(argv[1]);
    int num_ap = atoi(argv[2]);;
    double file_size = atoi(argv[3]);; // Mbyte

	int wifi_para = (argc-1)-3;

	if(double(wifi_para/num_ap)!=2){
		cerr << "Enter Parameters : <Bandwidth1> <bandwidth2> ... <LAMDA1> <LAMDA2>";
		return 1;
	}

    double lamda[num_ap]={0};
    double bandwidth[num_ap] = {0};

	for(int i=0;i<num_ap;i++)
	{
		bandwidth[i]=atoi(argv[4+i]);
		lamda[i]=atoi(argv[4+num_ap+i]);
	}

    double count1[num_ap] = {0};
    double count2[num_ap] = {0};
    double ti1[num_ap]={0};
    double ti2[num_ap]={0};
    cout << "*********** Select Best WiFi AP ***********" <<endl;
    cout << "************* Settings *************" << endl;
    cout << "Simulation Rounds :" << N_Simulation << endl;
    cout << "Bandwidth : " ;
    for(int i=0;i<num_ap;i++)
    {
        cout << "B["<<i<<"]=" <<bandwidth[i] << "MB/s ";
    }
    cout << endl;
    cout << "file_size = " << file_size << "Mbyte" <<endl;\

    cout << "lamda : " ;
    for(int i=0;i<num_ap;i++)
    {
         cout << " lamda["<<i<<"]=" <<lamda[i];
    }
    cout << endl;


    Random rm;
    double popsum  = 1;
    double arr[num_ap]={0};
    double mu[num_ap]={0};
    int ni[num_ap]={0};
    double Eni[num_ap]={0};
    double eeNi[num_ap]={0};
    // generate poisson random number...
    
    cout << "mu : " ;
    for(int i=0;i<num_ap;i++)
    {
         mu[i] = (double)(bandwidth[i]/file_size);
         cout << " mu["<<i<<"]=" <<mu[i];
    }
    cout << endl;

    cout << "generated ni : ";
    for(int i=0;i<num_ap;i++)
    {
        ni[i] = rm.poissonRandomNumber(lamda[i]);
        //ni[i] = rm.poissonRandomNumber(mu[i]);
        Eni[i]+=ni[i];
        cout <<" AP"<<i<<"="<<ni[i];
    }
    cout << endl;

    for(int n=0;n<N_Simulation;n++)
    {
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
        count1[index]+=1;

        //generate ti
        //double ti = rm.erlangRandomnumber(mu[index],ni[index]+1);
        //cout << "download time ti = " << ti << endl;

        //ti1[index]+=ti;

       //cout << "************** With WIFI Selection **************" <<endl; 
        //calculate min E[Ti]

        double eti[num_ap]={0};
        for(int i=0;i<num_ap;i++)
        {
            eti[i]= rm.erlangRandomnumber(mu[i],(ni[i]+1));
        }

        int selectwifi = FindMinIndex(eti,num_ap);

        count2[selectwifi]+=1;
        //double sti = rm.erlangRandomnumber(mu[selectwifi],ni[selectwifi]+1);
       // cout << "download time ti = " << sti << endl;
       // ti2[selectwifi]+=sti;
    }
    cout << "***********Simulation result****************"<<endl;
    cout << "==== without selection ====" <<endl;
    cout << "# of selected WiFi AP " <<endl;
    for(int i=0;i<num_ap;i++)
    {
        cout <<"AP["<<i<<"]="<<count1[i]<<endl;
    }

    cout << "************* Propability *************" << endl;
    for(int i=0;i<num_ap;i++)
    {
        double p = (double)(count1[i]/N_Simulation);
        double eti = (double)(ti1[i]/N_Simulation);
        cout <<"P[AP"<<i<<"]="<<p<<endl;
    }

    cout << "==== with selection ====" <<endl;
    cout << "# of selected WiFi AP " <<endl;
    for(int i=0;i<num_ap;i++)
    {
        cout <<"Select AP["<<i<<"]="<<count2[i]<<endl;
    }

    cout << "************* Propability *************" << endl;
    for(int i=0;i<num_ap;i++)
    {
        double p = (double)(count2[i]/N_Simulation);
        //double eti = (double)(ti2[i]/N_Simulation);
        cout <<"P[AP"<<i<<"]="<<p<<endl;
    }

    cout << "***********Math result****************"<<endl;

    ErlangDistribution es[num_ap]; // generate random variable

	for(int i=0;i<num_ap;i++){
		ErlangDistribution e (ni[i],mu[i]);
		es[i]= e;
	}

    double Mprob[num_ap]={0};
	double p2 = 0,psum =0;
    for(int i=0;i<num_ap;i++){
		// new function
       // cout << "shape = " <<es[i].getShape() << " lamda="<<es[i].getLamda()<<endl;
		p2 = getProbValue(es, i,1,num_ap);
        Mprob[i]=p2;
        cout << "P[AP"<<i<<" is the best]="<<p2<<endl;
        psum+=p2;
    }
    cout <<"Sum of P :"<< psum<<endl;

    return 0;
}
