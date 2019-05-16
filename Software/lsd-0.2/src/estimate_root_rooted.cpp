#include "estimate_root_rooted.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////estimate root without constraint (by using LD algorithm)///////////////////////////////////////////

double without_constraint_lambda(int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,double &lambda){
    //P: rooted tree
    //compute optimized solution without constraint (LD algorithm) with the position of root
	int r=Suc2[0];//r=r1
	int pr=Suc1[0];//pr=r2
	if (br==0) {
        lambda=0.5;
        B[r]=0;
        B[pr]=0;
        return without_constraint(n,m,P,Suc1,Suc2,B,V,T,rho,rho_min);
	}
	else{
        list<int> pos = postorder(P,Suc1,Suc2,n);
        double *W= new double[n];
        double *C = new double[n];
        double *X = new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho
        double *D= new double[n];//T[i]=D[i]+E[i]/rho
        double * E = new double[n];
        for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
            int i =  *iter;
            int s1 = Suc1[i];
            int s2 = Suc2[i];
            if (i==0){
                if (s1>=n && s2>=n){
                    D[i]=(T[s1]+T[s2])/2.;
                    E[i]=-br/2.;
                }
                else if (s1>=n){
                    D[i]=(T[s1]+D[s2])/2.;
                    E[i]=(E[s2]-br)/2.;
                }
                else if (s2>=n){
                    D[i]=(T[s2]+D[s1])/2.;
                    E[i]=(E[s1]-br)/2.;
                }
                else{
                    D[i]=(D[s1]+D[s2])/2.;
                    E[i]=(E[s1]+E[s2]-br)/2.;
                }
            }
            else if (i==r || i==pr){
                if (s1>=n && s2>=n){
                    D[i]=(T[s1]+T[s2])/2.;
                    E[i]=(-B[s1]-B[s2])/2.;
                }
                else if (s1>=n){
                    D[i]=(T[s1]+C[s2])/(2-W[s2]);
                    E[i]=(-B[s1]-B[s2]+X[s2])/(2-W[s2]);
                }
                else if (s2>=n){
                    D[i]=(T[s2]+C[s1])/(2-W[s1]);
                    E[i]=(-B[s1]-B[s2]+X[s1])/(2-W[s1]);
                }
                else{
                    D[i]=(C[s1]+C[s2])/(2-W[s1]-W[s2]);
                    E[i]=(-B[s1]-B[s2]+X[s1]+X[s2])/(2-W[s1]-W[s2]);
                }
            }
            else{
                if (s1>=n && s2>=n){
                    W[i]=1./3.;
                    C[i]=(T[s1]+T[s2])/3.;
                    X[i]=(-B[s1]-B[s2]+B[i])/3.;
                }
                else if (s1>=n){
                    W[i]=1./(3-W[s2]);
                    C[i]=(T[s1]+C[s2])/(3-W[s2]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s2])/(3-W[s2]);
                }
                else if (s2>=n){
                    W[i]=1./(3-W[s1]);
                    C[i]=(T[s2]+C[s1])/(3-W[s1]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s1])/(3-W[s1]);
                }
                else{
                    W[i]=1./(3-W[s1]-W[s2]);
                    C[i]=(C[s1]+C[s2])/(3-W[s1]-W[s2]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s1]+X[s2])/(3-W[s1]-W[s2]);
                }
            }
        }
        list<int> pre = preorder(P,Suc1,Suc2,n);
        for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
            int i = *iter;
            if (i!=r && i!=pr && i!=0){
                D[i]=W[i]*D[P[i]]+C[i];
                E[i]=W[i]*E[P[i]]+X[i];
            }
        }
        double M,N;//lambda=M*rho+N
        if (r>=n){
           M=(T[r]-D[pr])/(2*br);
           N=(-E[pr])/(2*br)+1./2.;
        }
        else if (pr>=n){
             M=(D[r]-T[pr])/(2*br);
             N=E[r]/(2*br)+1./2.;             
        }
        else{
            M=(D[r]-D[pr])/(2*br);
            N=(E[r]-E[pr])/(2*br)+1./2.;
        }
        double *G= new double[m+1];//
        double *H = new double[m+1];//
        if (r>=n){
           G[r]=M*br+D[0]-T[r];
           H[r]=N*br+E[0];
           G[pr]=-M*br+D[0]-D[pr];           
           H[pr]=(1-N)*br+E[0]-E[pr];
        }
        else if (pr>=n){
             G[r]=M*br+D[0]-D[r];
             H[r]=N*br+E[0]-E[r];
             G[pr]=-M*br+D[0]-T[pr];           
             H[pr]=(1-N)*br+E[0];                     
        }
        else{
             G[r]=M*br+D[0]-D[r];
             H[r]=N*br+E[0]-E[r];
             G[pr]=-M*br+D[0]-D[pr];           
             H[pr]=(1-N)*br+E[0]-E[pr];
        } 
        for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
            int i = *iter;
            if (i!=r && i!=pr && i!=0){
                G[i] = D[P[i]]-D[i];
                H[i] = B[i]+E[P[i]]-E[i];
            }
        }
        for (int i = n;i<=m;i++){
            if (i!=r && i!=pr){
                G[i] = D[P[i]]-T[i];
                H[i] = B[i]+E[P[i]];
            }
        }
        double a = 0;
        double b = 0;
        double c = 0;
        for (int i=1;i<=m;i++){
	  		a = a + G[i]*G[i];
	  		b = b + 2*G[i]*H[i];
	  		c = c + H[i]*H[i];
        }
        rho=-b/(2*a);
        if (rho<rho_min) rho=rho_min;
        lambda=M*rho+N;
        B[r]=lambda*br;
        B[pr]=(1-lambda)*br;
        for (int i=0;i<n;i++) T[i]=D[i]+E[i]/rho;
        delete[] W;
        delete[] C;
        delete[] X;
        delete[] D;
        delete[] E;
        delete[] G;
        delete[] H;
        if (lambda>1){
            lambda=1;
            B[r]=br;
            B[pr]=0;
            return without_constraint(n,m,P,Suc1,Suc2,B,V,T,rho,rho_min);
        }
        else if (lambda<0){
            lambda=0;
            B[r]=0;
            B[pr]=br;
            return without_constraint(n,m,P,Suc1,Suc2,B,V,T,rho,rho_min);
        }
        else{
            return phi(n,m,P,B,V,T,rho);
        }
    }
}

double estimate_root_without_constraint_local_rooted(int n,int m,int* &Suc1,int* &Suc2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda,double &rho,double rho_min){
    //P: rooted tree, recherche la nouvelle racine autour de l'ancien racine.
    //estimate root locally with LD algorithm for rooted tree////////////////////////
    cout<<"Re-estimating the root without constraints around the given root ... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    int s1=Suc1[0];
    int s2=Suc2[0];
    double br=reroot_rootedtree(m,s1,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    double ld,rhor;
    //cout<<"Optimizing the root position on the original branch "<<s1<<" ... ";
    cv[s1]=without_constraint_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
    //printf("%.10f\n",cv[s1]);
    cv[s2]=cv[s1];
    rho=rhor;
    lambda=ld;
    r=s1;
    phi1=cv[s1];
    list<int> next;
    if (s1<n){
        next.push_back(Suc1[s1]);
        next.push_back(Suc2[s1]);
    }
    if (s2<n){
        next.push_back(Suc1[s2]);
        next.push_back(Suc2[s2]);
    }
    while (!next.empty()){
        int i = next.back();
        br=reroot_rootedtree(m,i,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<i<<" ... ";
        cv[i]=without_constraint_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
        //printf("%.10f\n",cv[i]);
        if (cv[i]<cv[P[i]] && rhor>0){
            if (i<n){
                next.push_back(Suc1[i]);
                next.push_back(Suc2[i]);
            }
            if (cv[i]<phi1){
                phi1=cv[i];r=i;rho=rhor;lambda=ld;
            }
        }
        next.remove(i);
    }
    if (r==s1 || r==s2) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;    
    delete[] cv;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_without_constraint_rooted(int n,int m,int s1,int s2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda,double &rho,double rho_min){
    //P: rooted tree
    //estimate root with LD algorithm for rooted tree//////////////////////////////////////////
    cout<<"Re-estimating the root without constraints on all branches... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double ld,rhor;
    int y=1;
    double br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=without_constraint_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
    //printf("%.10f\n",phi1);
    lambda=ld;
    rho=rhor;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=without_constraint_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
            rho=rhor;
        }
        y++;
    }
    if (r==s1 || r==s2) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    return phi1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////estimate root by LD algorirthm (without constraint) in the case that the rate is known//////////////////////////////////
double without_constraint_lambda_rate(double rho,int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &lambda){
    //P: rooted tree
 	int r=Suc2[0];
	int pr=Suc1[0];
	if (br==0){
		lambda=0.5;
		B[r]=0;
		B[pr]=0;
		return without_constraint_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
	}
	else{
        list<int> pos = postorder(P,Suc1,Suc2,n);
        double *W= new double[n];
        double *C = new double[n];
        double *X = new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho
        double *D= new double[n];//T[i]=D[i]+E[i]/rho
        double * E = new double[n];
        for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
            int i =  *iter;
            int s1 = Suc1[i];
            int s2 = Suc2[i];
            if (i==0){
                if (s1>=n && s2>=n){
                    D[i]=(T[s1]+T[s2])/2.;
                    E[i]=-br/2.;
                }
                else if (s1>=n){
                    D[i]=(T[s1]+D[s2])/2.;
                    E[i]=(E[s2]-br)/2.;
                }
                else if (s2>=n){
                    D[i]=(T[s2]+D[s1])/2.;
                    E[i]=(E[s1]-br)/2.;
                }
                else{
                    D[i]=(D[s1]+D[s2])/2.;
                    E[i]=(E[s1]+E[s2]-br)/2.;
                }
            }
            else if (i==r || i==pr){
                if (s1>=n && s2>=n){
                    D[i]=(T[s1]+T[s2])/2.;
                    E[i]=(-B[s1]-B[s2])/2.;
                }
                else if (s1>=n){
                    D[i]=(T[s1]+C[s2])/(2-W[s2]);
                    E[i]=(-B[s1]-B[s2]+X[s2])/(2-W[s2]);
                }
                else if (s2>=n){
                    D[i]=(T[s2]+C[s1])/(2-W[s1]);
                    E[i]=(-B[s1]-B[s2]+X[s1])/(2-W[s1]);
                }
                else{
                    D[i]=(C[s1]+C[s2])/(2-W[s1]-W[s2]);
                    E[i]=(-B[s1]-B[s2]+X[s1]+X[s2])/(2-W[s1]-W[s2]);
                }
            }
            else{
                if (s1>=n && s2>=n){
                    W[i]=1./3.;
                    C[i]=(T[s1]+T[s2])/3.;
                    X[i]=(-B[s1]-B[s2]+B[i])/3.;
                }
                else if (s1>=n){
                    W[i]=1./(3-W[s2]);
                    C[i]=(T[s1]+C[s2])/(3-W[s2]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s2])/(3-W[s2]);
                }
                else if (s2>=n){
                    W[i]=1./(3-W[s1]);
                    C[i]=(T[s2]+C[s1])/(3-W[s1]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s1])/(3-W[s1]);
                }
                else{
                    W[i]=1./(3-W[s1]-W[s2]);
                    C[i]=(C[s1]+C[s2])/(3-W[s1]-W[s2]);
                    X[i]=(-B[s1]-B[s2]+B[i]+X[s1]+X[s2])/(3-W[s1]-W[s2]);
                }
            }
        }
        list<int> pre = preorder(P,Suc1,Suc2,n);
        for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
            int i = *iter;
            if (i!=r && i!=pr && i!=0){
                D[i]=W[i]*D[P[i]]+C[i];
                E[i]=W[i]*E[P[i]]+X[i];
            }
        }
        double M,N;//lambda=M*rho+N
        if (r>=n){
           M=(T[r]-D[pr])/(2*br);
           N=(-E[pr])/(2*br)+1./2.;
        }
        else if (pr>=n){
             M=(D[r]-T[pr])/(2*br);
             N=E[r]/(2*br)+1./2.;             
        }
        else{
            M=(D[r]-D[pr])/(2*br);
            N=(E[r]-E[pr])/(2*br)+1./2.;
        }
        lambda=M*rho+N;
        B[r]=lambda*br;
        B[pr]=(1-lambda)*br;
        for (int i=0;i<n;i++) T[i]=D[i]+E[i]/rho;
        delete[] W;
        delete[] C;
        delete[] X;        
        delete[] D;
        delete[] E;
        if (lambda>1){
            lambda=1;
            B[r]=br;
            B[pr]=0;
            return without_constraint_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
        }
        else if (lambda<0){
            lambda=0;
            B[r]=0;
            B[pr]=br;
            return without_constraint_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
        }
        else{
            return phi(n,m,P,B,V,T,rho);
        }
    }
}

double estimate_root_without_constraint_local_rate_rooted(double rho,int n,int m,int* &Suc1,int* &Suc2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda){
    //P: rooted tree, recherche la nouvelle racine autour de l'ancien racine.
    //estimate locally the root without constraint, the rate is given
    cout<<"Re-estimating the root without constraints around the given root ... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    int s1=Suc1[0];
    int s2=Suc2[0];
    double br=reroot_rootedtree(m,s1,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    double ld;
    //cout<<"Optimizing the root position on the original branch "<<s1<<" ... ";
    cv[s1]=without_constraint_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
    //printf("%.10f\n",cv[s1]);
    cv[s2]=cv[s1];
    lambda=ld;
    r=s1;
    phi1=cv[s1];
    list<int> next;
    if (s1<n){
        next.push_back(Suc1[s1]);
        next.push_back(Suc2[s1]);
    }
    if (s2<n){
        next.push_back(Suc1[s2]);
        next.push_back(Suc2[s2]);
    }
    while (!next.empty()){
        int i = next.back();
        br=reroot_rootedtree(m,i,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<i<<" ... ";
        cv[i]=without_constraint_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
        //printf("%.10f\n",cv[i]);
        if (cv[i]<cv[P[i]]){
            if (i<n){
                next.push_back(Suc1[i]);
                next.push_back(Suc2[i]);
            }
            if (cv[i]<phi1){
                phi1=cv[i];r=i;lambda=ld;
            }
        }
        next.remove(i);
    }
    if (r==Suc1[0] || r==Suc2[0]) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] cv;
    delete[] P_new;
    delete[] Suc1;
    delete[] Suc2;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_without_constraint_rate_rooted(double rho,int n,int m,int s1,int s2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda){
    //P: rooted tree
    cout<<"Re-estimating the root without constraints on all branches ... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double ld;
    int y=1;
    double br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=without_constraint_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
    //printf("%.10f\n",phi1);
    lambda=ld;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=without_constraint_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
        }
        y++;
    }
    if (r==s1 || r==s2) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////estimate root with constraint (by using QPD algorithm)/////////////////////////////////////////
double starting_point_lambda(int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,double &lambda,list<int>  &active_set){
    //compute the starting point//////////////
    int s1 = Suc1[0];
    int s2 = Suc2[0];
    without_constraint_lambda(n,m,B[s1]+B[s2],P,Suc1,Suc2,B,V,T,rho,rho_min,lambda);
    list<int> pos = postorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (T[i]>T[s1] || T[i]>T[s2]){
            if (T[s1]<T[s2]){
                T[i]=T[s1];
                active_set.push_back(s1);
            }
            else{
                T[i]=T[s2];
                active_set.push_back(s2);
            }
        }
    }
    return phi(n,m,P,B,V,T,rho);
}


bool conditions_lambda(int m,list<double>& ldLagrange,int* & P,double* & T,double lambda){
    //stop condition
    for (list<double>::iterator iter = ldLagrange.begin();iter!=ldLagrange.end();iter++){
        if (*iter < 0){
            return false;
        }
    }
    for (int i=1;i<=m;i++) {
        if (T[i]<T[P[i]]-(10e-15)) return false;
    }
    if (lambda<0) return false;
    if (lambda>1) return false;
    return true;
}

list<double> without_constraint_active_set_lambda(int n,int m,double br,int* &P,bool* & flag,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,double &lambda,list<int> &active_set,double &phi1){
    //P: rooted tree with unknown lambda, br: length of the branch where stiuated the root
    //this methods computes the optimized solution of the Lagrange function
    int r = Suc2[0];
    int pr = Suc1[0];
    if (br==0) {
        lambda=0.5;
        B[r]=0;
        B[pr]=0;
        return without_constraint_active_set(n,m,flag,P,Suc1,Suc2,B,V,T,rho,rho_min,active_set,phi1);
    }
    else{
        double* ldLagrange = new double[active_set.size()];
        list<double> ld;
        list<int>* SucL = new list<int>[n];
        list<int>* SucI = new list<int>[n];
        
        int *Pre= new int[m+1];
        for (int i=0;i<=m;i++) Pre[i]=-1;
        list<int> ls;
        for (int i=0;i<=m+2;i++) flag[i]=true;//flag[i]=false means i is in the active set
        bool *leaf= new bool[m+1];
        for (int i=0;i<n;i++) leaf[i]=false;
        for (int i=n;i<=m;i++) leaf[i]=true;
        //flag
        for (list<int>::iterator iter=active_set.begin();iter!=active_set.end();iter++){
            int i = *iter;
            flag[i]=false;
            if (i>=n && i!=m+1 && i!=m+2) ls.push_back(i);//m+1: lambda>=0, m+2: lambda<=1
        }
        
        stack<int>* feuilles = new stack<int>[ls.size()];
        
        list<int> top;
        //leaf
        computeFeuilles(n,m,P,Suc1,Suc2,T,flag,leaf,feuilles,active_set);
        for (int i=0;i<n;i++){
            int s1= Suc1[i];
            int s2= Suc2[i];
            if (flag[i] && !leaf[i] && (!flag[s1] || !flag[s2]))
                top.push_back(i);
        }
        
        list<int>* internal = new list<int>[top.size()];
        
        reduceTree(n,m,P,Pre,Suc1,Suc2,SucL,SucI,flag,leaf,internal,active_set);
        
        list<int> pos = postorder(P,Suc1,Suc2,n);
        list<int> pre = preorder(P,Suc1,Suc2,n);
        
        
        if (!flag[m+1]){
            B[r]=0;
            B[pr]=br;
            lambda=0;
            return without_constraint_active_set(n,m,flag,P,Suc1,Suc2,B,V,T,rho,rho_min,active_set,phi1);
        }
        else if (!flag[m+2]){
            B[r]=br;
            B[pr]=0;
            lambda=1;
            return without_constraint_active_set(n,m,flag,P,Suc1,Suc2,B,V,T,rho,rho_min,active_set,phi1);
        }
        else{
            double *W= new double[n];
            double *C = new double[n];
            double *X = new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho
            double *D= new double[n];//T[i]=D[i]+E[i]/rho
            double * E = new double[n];
            for (int i=0;i<n;i++){
                W[i]=0;C[i]=0;X[i]=0;D[i]=0;E[i]=0;
            }
            for (list<int>::iterator it = pos.begin();it!=pos.end();it++){
				int i = *it;
				if (flag[i] && !leaf[i]){
					int p = Pre[i];
					if (p==-1){
						if (flag[r] && flag[pr]){
                            D[i]=0;
                            E[i]=-br/2.;
                            for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                                D[i]+=T[*iter]/2.;
                            }
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                                D[i]+=D[*iter]/2.;
                                E[i]+=E[*iter]/2.;
                            }
						}
						else if (flag[r] || flag[pr]){
							double coefs=1./2.;
							double ctemp=0;
							double xtemp=-br/2.;
							double k=r;
							if (flag[pr]) k=pr;
							for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
								if (*iter!=k){
									coefs+=1;
									ctemp+=T[*iter];
									xtemp-=B[*iter];
								}
								else{
                                    ctemp+=T[*iter]/2.;
                                }
							}
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
								if (*iter!=k){
                                    coefs+=(1-W[*iter]);
                                    ctemp+=C[*iter];
                                    xtemp+=(X[*iter]-B[*iter]);
								}
								else{
									coefs+=-W[*iter]/2.;
									ctemp+=C[*iter]/2.;
									xtemp+=X[*iter]/2.;
                                }
							}
							D[i]=ctemp/coefs;
							E[i]=xtemp/coefs;
						}
						else{//!flag[r] && !flag[pr]
							double coefs=0;
							double xtemp=0;
							double ctemp=0;
							for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
								coefs+=1;
								ctemp+=T[*iter];
								xtemp-=B[*iter];
							}
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
								coefs+=(1-W[*iter]);
								ctemp+=C[*iter];
								xtemp+=(X[*iter]-B[*iter]);
							}
							D[i]=ctemp/coefs;
							E[i]=xtemp/coefs;
						}
					}
					else if ((i==r || i==pr) && flag[r] && flag[pr]){
						double coefs=0;
						double xtemp=0;
						double ctemp=0;
						for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
							coefs+=1;
							ctemp+=T[*iter];
							xtemp-=B[*iter];
						}
						for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
							coefs+=(1-W[*iter]);
							ctemp+=C[*iter];
							xtemp+=(X[*iter]-B[*iter]);
						}
						D[i]=ctemp/coefs;
						E[i]=xtemp/coefs;
					}
					else if ((i==r && !flag[pr])||(i==pr && !flag[r])){
						double coefs=1./2.;
						double wtemp=1./2.;
						double ctemp=0;
						double xtemp=br/2.;
						for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
							coefs+=1;
							ctemp+=T[*iter];
							xtemp-=B[*iter];
					    }
						for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
							coefs+=(1-W[*iter]);
							ctemp+=C[*iter];
							xtemp+=(X[*iter]-B[*iter]);
                        }
                        if (!leaf[p]){
                            W[i]=wtemp/coefs;
                            C[i]=ctemp/coefs;
                            X[i]=xtemp/coefs;
                        }
                        else{
                            ctemp+=T[p]/2.;
                            D[i]=ctemp/coefs;
                            E[i]=xtemp/coefs;
                        }
					}
					else{
                        double coefs=1;
                        double wtemp=1;
                        double ctemp=0;
                        double xtemp=B[i];
                        for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                            coefs+=1;
                            ctemp+=T[*iter];
                            xtemp-=B[*iter];
                        }
                        for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                            coefs+=(1-W[*iter]);
                            ctemp+=C[*iter];
                            xtemp+=(X[*iter]-B[*iter]);
                        }
                        if (!leaf[p]){
                            W[i]=wtemp/coefs;
							C[i]=ctemp/coefs;
							X[i]=xtemp/coefs;
                        }
                        else{
                            ctemp+=T[p];
                            D[i]=ctemp/coefs;
                            E[i]=xtemp/coefs;
                        }
					}
				}
			}
            for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
                int i = *iter;
                if ((Pre[i]!=-1 && flag[i] && !leaf[i] && !leaf[Pre[i]]) && !((i==r || i==pr) && flag[r] && flag[pr])){
                    D[i]=W[i]*D[Pre[i]]+C[i];
                    E[i]=W[i]*E[Pre[i]]+X[i];
                }
            }
            for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
                int i = *iter;
                int s1 = Suc1[i];
                int s2 = Suc2[i];
                if (!flag[s1] && !leaf[s1]) {
                    D[s1]=D[i];
                    E[s1]=E[i];
                }
                if (!flag[s2] && !leaf[s2]) {
                    D[s2]=D[i];
                    E[s2]=E[i];
                }
            }
            //lambda = rho*Dlambda+Elambda
            double Dlambda=0;
            double Elambda=1./2.;
            if (flag[r] && flag[pr]){
                if (leaf[r]) Dlambda+=T[r]/(2*br);
                else {
                    Dlambda+=D[r]/(2*br);
                    Elambda+=E[r]/(2*br);
                }
                if (leaf[pr]) Dlambda-=T[pr]/(2*br);
                else {
                    Dlambda-=D[pr]/(2*br);
                    Elambda-=E[pr]/(2*br);
                }
            }
            else if (flag[r] && !flag[pr]){
                if (leaf[r]) Dlambda+=T[r]/(2*br);
                else {
                    Dlambda+=D[r]/(2*br);
                    Elambda+=E[r]/(2*br);
                }
                if (leaf[0]) Dlambda-=T[0]/(2*br);
                else {
                    Dlambda-=D[0]/(2*br);
                    Elambda-=E[0]/(2*br);
                }
            }
            else if (!flag[r] && flag[pr]){
                if (leaf[0]) Dlambda+=T[0]/(2*br);
                else {
                    Dlambda+=D[0]/(2*br);
                    Elambda+=E[0]/(2*br);
                }
                if (leaf[pr]) Dlambda-=T[pr]/(2*br);
                else {
                    Dlambda-=D[pr]/(2*br);
                    Elambda-=E[pr]/(2*br);
                }
            }
            double *F= new double[m+1];
            double *G= new double[m+1];//T[P[i]]-T[i]=F[i]+G[i]/rho;
            for (int i=0;i<=m;i++){
                F[i]=0;G[i]=0;
            }
            
            for (int i=1;i<=m;i++){
                if (leaf[P[i]]) F[i]+=T[P[i]];
                else {
                    F[i]+=D[P[i]];
                    G[i]+=E[P[i]];
                }
                if (leaf[i]) F[i]-=T[i];
                else {
                    F[i]-=D[i];
                    G[i]-=E[i];
                }
            }
            double a = 0;
            double b = 0;
            double c = 0; //a*rho*rho+b*rho*c=0
            for (int i=1;i<=m;i++){
                if (i!=r && i!=pr){
                    a+=F[i]*F[i];
                    b+=2*F[i]*G[i]+B[i]*F[i];
                    c+=G[i]*G[i]+B[i]*G[i];
                }
            }
            a+=F[r]*F[r]+F[pr]*F[pr]+br*Dlambda*(F[r]-F[pr]);
            b+=2*F[r]*G[r]+2*F[pr]*G[pr]+br*Dlambda*(G[r]-G[pr])+Elambda*br*(F[r]-F[pr])+br*F[pr];
            c+=G[r]*G[r]+G[pr]*G[pr]+br*Elambda*(G[r]-G[pr])+br*G[pr];
            double delta = b*b-4*a*c;
            if (delta>=0){
                rho = (-b+sqrt(delta))/(2*a);
                if (rho<rho_min) rho=rho_min;
            }
            else {
                rho=rho_min;
            }
            lambda=Dlambda*rho+Elambda;
            for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
                int i = *iter;
                if (flag[i] && !leaf[i]) T[i] = D[i]+E[i]/rho;
                int s1 = Suc1[i];
                int s2 = Suc2[i];
                if (!flag[s1]) T[s1]=T[i];
                if (!flag[s2]) T[s2]=T[i];
            }
            
            B[r]=br*lambda;
            B[pr]=br-br*lambda;
            phi1 = phi(n,m,P,B,V,T,rho);
            delete[] W;
            delete[] C;
            delete[] X;
            delete[] D;
            delete[] E;
            delete[] F;
            delete[] G;
        }
        
        int *as= new int[m+3];
        for (int i=0;i<=m+2;i++) as[i]=-1;
        int count = 0;
        for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
            as[*iter]=count;
            count++;
        }
        bool *bl= new bool[m+1];//bl[i]=true means ldLagrange[i] has been computed
        for (int i=0;i<=m;i++) bl[i]=false;
        for (int y=0;y<ls.size();y++){
            while (!feuilles[y].empty()){
                int i=feuilles[y].top();
                if (!flag[i]){
                    feuilles[y].pop();
                    int p=P[i];
                    int s1=Suc1[i];
                    int s2=Suc2[i];
                    if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                            ldLagrange[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                        if (flag[s1] && flag[s2]) bl[i]=true;
                        if (!flag[s1] && bl[s1] && flag[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]];bl[i]=true;}
                        if (!flag[s2] && bl[s2] && flag[s1]) {ldLagrange[as[i]]+=ldLagrange[as[s2]];bl[i]=true;}
                        if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]]+ldLagrange[as[s2]];bl[i]=true;}
                    }
                    else{
                        int j=Suc1[p];
                        if (j==i)
                            j=Suc2[p];
                        ldLagrange[as[i]]=-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i]-2*rho*(B[j]-rho*T[j]+rho*T[p])/V[j];
                        if (P[p]!=-1)
                            ldLagrange[as[i]]+=2*rho*(B[p]-rho*T[p]+rho*T[P[p]])/V[p];
                        if (flag[j] && flag[p]) bl[i]=true;
                        if (flag[p] && !flag[j] && bl[j]) {ldLagrange[as[i]]-=ldLagrange[as[j]];bl[i]=true;}
                        if (!flag[p] && bl[p] && flag[j]) {ldLagrange[as[i]]+=ldLagrange[as[p]];bl[i]=true;}
                        if (!flag[j] && bl[j] && !flag[p] && bl[p]){ldLagrange[as[i]]+=ldLagrange[as[p]]-ldLagrange[as[j]];bl[i]=true;}
                    }
                }
            }
        }
        for (int y=0;y<top.size();y++){
            for (list<int>::iterator iter = internal[y].begin();iter!=internal[y].end();iter++){
                int i= *iter;
                int p = P[i];
                int s1=Suc1[i];
                int s2=Suc2[i];
                if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                        ldLagrange[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                    if (flag[s1] && flag[s2]) bl[i]=true;
                    if (!flag[s1] && bl[s1] && flag[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]];bl[i]=true;}
                    if (!flag[s2] && bl[s2] && flag[s1]) {ldLagrange[as[i]]+=ldLagrange[as[s2]];bl[i]=true;}
                    if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]]+ldLagrange[as[s2]];bl[i]=true;}
                }
            }
        }
        if (!flag[m+1]){
            ldLagrange[as[m+1]]=2*br*(br*lambda-rho*T[r]+rho*T[0])-2*br*(br-br*lambda-rho*T[pr]+rho*T[0]);
        }
        if (!flag[m+2]){
            ldLagrange[as[m+2]]=-2*br*(br*lambda-rho*T[r]+rho*T[0])+2*br*(br-br*lambda-rho*T[pr]+rho*T[0]);
        }
        
        for (int i=0;i<active_set.size();i++){
            if (myabs(ldLagrange[i])<(10e-15))
                ldLagrange[i]=0;
            ld.push_back(ldLagrange[i]);
        }
        delete[] Pre;
        delete[]  leaf;
        delete[] SucL;
        delete[] SucI;
        delete[] feuilles;
        delete[] internal;
        delete[] as;
        delete[] bl;
        delete[] ldLagrange;
        return ld;
    }
}

double with_constraint_active_set_lambda(int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,double &lambda){
    //compute the optimized solution with constraint and with the position of root
    list<int> active_set;
    double phi1;
    phi1 = starting_point_lambda(n,m,br,P,Suc1,Suc2,B,V,T,rho,rho_min,lambda,active_set);
    double *T_new= new double[m+1];
    double *dir= new double[n];
    double alpha;
    for (int i=n;i<=m;i++) T_new[i]=T[i];
    bool *flag= new bool[m+3];
    double lambda_new;
    list<double> ldLagrange = without_constraint_active_set_lambda(n,m,br,P,flag,Suc1,Suc2,B,V,T_new,rho,rho_min,lambda_new,active_set,phi1);
    while (!conditions_lambda(m,ldLagrange,P,T_new,lambda_new)){
        for (int i=0;i<n;i++)  dir[i]=T_new[i]-T[i];
        double dirlambda = lambda_new-lambda;
        alpha=1;
        int as=0;
        double a;
        if (dirlambda<0){
            a = lambda/(-dirlambda);
            if (a<alpha){
                alpha=a;
                as = m+1;
            }
        }
        if (dirlambda>0){
            a = (1-lambda)/dirlambda;
            if (a<alpha){
                alpha=a;
                as=m+2;
            }
        }
        for (int i=n;i<=m;i++){
            if (flag[i]){
                int p = P[i];
                if (dir[p]>0){
                    a = (T[i]-T[p])/dir[p];
                    if (a<alpha){
                        alpha = a;
                        as = i;
                    }
                }
            }
        }
        for (int i=0;i<n;i++){
            if (flag[i]){
                int p = P[i];
                if (p!=-1){
                    if (dir[p]>dir[i]){
                        a = (T[i]-T[p])/(dir[p]-dir[i]);
                        if (a<alpha){
                            alpha = a;
                            as = i;
                        }
                    }
                }
            }
        }
        for (int i=0;i<n;i++) T[i]=T[i]+alpha*dir[i];
        lambda=lambda+alpha*dirlambda;
        int asrm= remove_ne_lambda(ldLagrange,active_set);
        if (asrm!=-1){
            active_set.remove(asrm);
        }
        if (as!=0){
            active_set.push_back(as);
        }
        ldLagrange = without_constraint_active_set_lambda(n,m,br,P,flag,Suc1,Suc2,B,V,T_new,rho,rho_min,lambda_new,active_set,phi1);
    }
    for (int i=0;i<n;i++) T[i]=T_new[i];
    lambda=lambda_new;
    delete[] T_new;
    delete[] dir;
    delete[] flag;
    return phi1;
}

double estimate_root_with_constraint_local_rooted(int n,int m,int* &Suc1,int* &Suc2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda,double &rho,double rho_min){
    //P: rooted tree, recherche la nouvelle racine autour de l'ancienne racine.
    /////////////estimate root locally with QPD algorithm for rooted tree////////////////////////////////////////////////////////////    
    cout<<"Re-estimating the root with constraints around the given root ... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    int s1=Suc1[0];
    int s2=Suc2[0];
    double br=reroot_rootedtree(m,s1,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    double ld,rhor;
    //cout<<"Optimizing the root position on the original branch "<<s1<<" ... ";
    cv[s1]=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
    //printf("%.10f\n",cv[s1]);
    cv[s2]=cv[s1];
    rho=rhor;
    lambda=ld;
    r=s1;
    phi1=cv[s1];
    list<int> next;
    if (s1<n){
        next.push_back(Suc1[s1]);
        next.push_back(Suc2[s1]);
    }
    if (s2<n){
        next.push_back(Suc1[s2]);
        next.push_back(Suc2[s2]);
    }
    while (!next.empty()){
        int i = next.back();
        br=reroot_rootedtree(m,i,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<i<<" ... ";
        cv[i]=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
        //printf("%.10f\n",cv[i]);
        if (cv[i]<cv[P[i]]){
            if (i<n){
                next.push_back(Suc1[i]);
                next.push_back(Suc2[i]);
            }
            if (cv[i]<phi1){
                phi1=cv[i];r=i;rho=rhor;lambda=ld;
            }
        }
        next.remove(i);
    }
    if (r==Suc1[0] || r==Suc2[0]) cout<<"The new root is on the original branch."<<endl; 
    else cout<<"The new root is on the branch "<<r<<endl;   
    delete[] cv;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_with_constraint_fast_rooted(int n,int m,int s01,int s02,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda,double &rho,double rho_min){
    //P: rooted tree, oublier la racine, recherche la nouvelle racine sur toutes les branches
    /////////////estimate root on all branches with QPD algorithm for rooted tree////////////////////////////////////////////////////////////
    cout<<"Re-estimating the root with constraints on all branches using fast method."<<endl;
    double ld,rhor;
    double phi1=estimate_root_without_constraint_rooted(n,m,s01,s02,P,B,V,T,r,ld,rhor,rho_min);
    cout<<"Re-estimating the root with constraints around the given root ... ";
    //cout<<"Optimizing the root position on the branch "<<r<<" ... ";   
    int* P_new = new int[m+1];
    int* Suc1_new= new int[n];
    int* Suc2_new= new int[n];  
    double* B_new= new double[m+1];
    double* T_new= new double[m+1];  
    int* P_ref = new int[m+1];
    int* tab = new int[m+1];;
    double br=reroot_rootedtree(m,r,P,s01,s02,B,P_new,B_new,P_ref,tab);
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    computeSuc(P_new, Suc1_new, Suc2_new, m+1, n);  
    cv[r]=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
    //printf("%.10f\n",cv[r]);
    list<int> next;
    phi1=cv[r];
    lambda=ld;
    rho=rhor;
    int* Suc1_ref = new int[n];
    int* Suc2_ref = new int[n];
    computeSuc(P_ref,Suc1_ref,Suc2_ref,m+1,n);
    int s1=Suc1_ref[0];
    int s2=Suc2_ref[0];
    if (Suc1_ref[0]<n){
        next.push_back(Suc1_ref[s1]);
        next.push_back(Suc2_ref[s1]);
    }
    if (Suc2_ref[0]<n){
        next.push_back(Suc1_ref[s2]);
        next.push_back(Suc2_ref[s2]);
    }
    while (!next.empty()){
        int i = next.back();
        int e = tab[i];
        br=reroot_rootedtree(m,e,P,s01,s02,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<e<<" ... ";
        cv[e]=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
        //printf("%.10f\n",cv[e]);
        if (cv[e]<cv[tab[P_ref[i]]]){
            if (i<n){
                next.push_back(Suc1_ref[i]);
                next.push_back(Suc2_ref[i]);
            }
            if (cv[i]<phi1){
                phi1=cv[e];r=e;rho=rhor;lambda=ld;
            }
        }
        next.remove(i);
    }
    if (r==s01 || r==s02) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] cv;
    delete[] P_new;
    delete[] P_ref;
    delete[] tab;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] Suc1_ref;
    delete[] Suc2_ref;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}


double estimate_root_with_constraint_rooted(int n,int m,int s1,int s2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda,double &rho,double rho_min){
    //P: rooted tree
    //estimate root with QPD algorithm for rooted tree//////////////////////////////////
    cout<<"Re-estimating the root using constrained mode on all branches ... ";
    double phi1=-1;
    int y=1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double ld,rhor;
    double br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
    //printf("%.10f\n",phi1);
    lambda=ld;
    rho=rhor;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=with_constraint_active_set_lambda(n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,rhor,rho_min,ld);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
            rho=rhor;
        }
        y++;
    }
    if (r==s1 || r==s2) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////estimate root by QPD algorithm (with constraint) in the case that the rate is known///////////////////////////
double starting_point_lambda_rate(double rho,int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &lambda,list<int> & active_set){
    int s1 = Suc1[0];
    int s2 = Suc2[0];
    without_constraint_lambda_rate(rho,n,m,B[s1]+B[s2],P,Suc1,Suc2,B,V,T,lambda);
    list<int> pos = postorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (T[i]>T[s1] || T[i]>T[s2]){
            if (T[s1]<T[s2]){
                T[i]=T[s1];
                active_set.push_back(s1);
            }
            else{
                T[i]=T[s2];
                active_set.push_back(s2);
            }
        }
    }
    return phi(n,m,P,B,V,T,rho);
}

list<double> without_constraint_active_set_lambda_rate(double rho,int n,int m,double br,int* & P,bool* & flag,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &lambda,list<int>& active_set,double &phi1){
    //P: rooted tree with unknown lambda, br: length of the branch where stiuated the root
	int r = Suc2[0];
	int pr = Suc1[0];
	if (br==0) {
		lambda=0.5;
		B[r]=0;
		B[pr]=0;
		return without_constraint_active_set_rate(rho,n,m,flag,P,Suc1,Suc2,B,V,T,active_set,phi1);
	}
	else{
        double* ldLagrange = new double[active_set.size()];
        list<double> ld;
        
        list<int>* SucL = new list<int>[n];
        list<int>* SucI = new list<int>[n];
        
        int *Pre= new int[m+1];
        for (int i=0;i<=m;i++) Pre[i]=-1;
        list<int> ls;
        for (int i=0;i<=m+2;i++) flag[i]=true;
        for (int i=0;i<=m;i++) flag[i]=true;
        bool *leaf= new bool[m+1];
        for (int i=0;i<n;i++) leaf[i]=false;
        for (int i=n;i<=m;i++) leaf[i]=true;
        //flag
        for (list<int>::iterator iter=active_set.begin();iter!=active_set.end();iter++){
            int i = *iter;
            flag[i]=false;
            if (i>=n && i!=m+1 && i!=m+2) ls.push_back(i);//m+1: lambda>=0, m+2: lambda<=1
        }
        
        // stack<int> feuilles[ls.size()];
        stack<int>* feuilles = new stack<int>[ls.size()];
        
        list<int> top;
        //leaf
        computeFeuilles(n,m,P,Suc1,Suc2,T,flag,leaf,feuilles,active_set);
        for (int i=0;i<n;i++){
            int s1= Suc1[i];
            int s2= Suc2[i];
            if (flag[i] && !leaf[i] && (!flag[s1] || !flag[s2]))
                top.push_back(i);
        }
        
        // list<int> internal[top.size()];
        list<int>* internal = new list<int>[top.size()];
        
        reduceTree(n,m,P,Pre,Suc1,Suc2,SucL,SucI,flag,leaf,internal,active_set);
        
        list<int> pos = postorder(P,Suc1,Suc2,n);
        list<int> pre = preorder(P,Suc1,Suc2,n);
        
        
        if (!flag[m+1]){
            B[r]=0;
            B[pr]=br;
            lambda=0;
            return without_constraint_active_set_rate(rho,n,m,flag,P,Suc1,Suc2,B,V,T,active_set,phi1);
        }
        else if (!flag[m+2]){
            B[r]=br;
            B[pr]=0;
            lambda=1;
            return without_constraint_active_set_rate(rho,n,m,flag,P,Suc1,Suc2,B,V,T,active_set,phi1);
        }
        else{
            double *W= new double[n];
            double *C = new double[n];
            double *X = new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho
            double *D= new double[n];//T[i]=D[i]+E[i]/rho
            double * E = new double[n];
            for (int i=0;i<n;i++){
                W[i]=0;C[i]=0;X[i]=0;D[i]=0;E[i]=0;
            }
            for (list<int>::iterator it = pos.begin();it!=pos.end();it++){
				int i = *it;
				if (flag[i] && !leaf[i]){
					int p = Pre[i];
					if (p==-1){
						if (flag[r] && flag[pr]){
                            D[i]=0;
                            E[i]=-br/2.;
                            for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                                D[i]+=T[*iter]/2;
                            }
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                                D[i]+=D[*iter]/2;
                                E[i]+=E[*iter]/2;
                            }
						}
						else if (flag[r] || flag[pr]){
							double coefs=1./2;
							double ctemp=0;
							double xtemp=-br/2.;
							double k=r;
							if (flag[pr]) k=pr;
							for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
								if (*iter!=k){
									coefs+=1;
									ctemp+=T[*iter];
									xtemp-=B[*iter];
								}
								else{
                                    ctemp+=T[*iter]/2.;
                                }
							}
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
								if (*iter!=k){
                                    coefs+=(1-W[*iter]);
                                    ctemp+=C[*iter];
                                    xtemp+=(X[*iter]-B[*iter]);
								}
								else{
									coefs+=-W[*iter]/2.;
									ctemp+=C[*iter]/2.;
									xtemp+=X[*iter]/2.;
                                }
							}
							D[i]=ctemp/coefs;
							E[i]=xtemp/coefs;
						}
						else{//!flag[r] && !flag[pr]
							double coefs=0;
							double xtemp=0;
							double ctemp=0;
							for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
								coefs+=1;
								ctemp+=T[*iter];
								xtemp-=B[*iter];
							}
							for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
								coefs+=(1-W[*iter]);
								ctemp+=C[*iter];
								xtemp+=(X[*iter]-B[*iter]);
							}
							D[i]=ctemp/coefs;
							E[i]=xtemp/coefs;
						}
					}
					else if ((i==r || i==pr) && flag[r] && flag[pr]){
						double coefs=0;
						double xtemp=0;
						double ctemp=0;
						for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
							coefs+=1;
							ctemp+=T[*iter];
							xtemp-=B[*iter];
						}
						for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
							coefs+=(1-W[*iter]);
							ctemp+=C[*iter];
							xtemp+=(X[*iter]-B[*iter]);
						}
						D[i]=ctemp/coefs;
						E[i]=xtemp/coefs;
					}
					else if ((i==r && !flag[pr])||(i==pr && !flag[r])){
						double coefs=1./2.;
						double wtemp=1./2.;
						double ctemp=0;
						double xtemp=br/2.;
						for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
							coefs+=1;
							ctemp+=T[*iter];
							xtemp-=B[*iter];
					    }
						for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
							coefs+=(1-W[*iter]);
							ctemp+=C[*iter];
							xtemp+=(X[*iter]-B[*iter]);
                        }
                        if (!leaf[p]){
                            W[i]=wtemp/coefs;
                            C[i]=ctemp/coefs;
                            X[i]=xtemp/coefs;
                        }
                        else{
                            ctemp+=T[p]/2.;
                            D[i]=ctemp/coefs;
                            E[i]=xtemp/coefs;
                        }
					}
					else{
                        double coefs=1;
                        double wtemp=1;
                        double ctemp=0;
                        double xtemp=B[i];
                        for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                            coefs+=1;
                            ctemp+=T[*iter];
                            xtemp-=B[*iter];
                        }
                        for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                            coefs+=(1-W[*iter]);
                            ctemp+=C[*iter];
                            xtemp+=(X[*iter]-B[*iter]);
                        }
                        if (!leaf[p]){
                            W[i]=wtemp/coefs;
							C[i]=ctemp/coefs;
							X[i]=xtemp/coefs;
                        }
                        else{
                            ctemp+=T[p];
                            D[i]=ctemp/coefs;
                            E[i]=xtemp/coefs;
                        }
					}
				}
			}
            for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
                int i = *iter;
                if ((Pre[i]!=-1 && flag[i] && !leaf[i] && !leaf[Pre[i]]) && !((i==r || i==pr) && flag[r] && flag[pr])){
                    D[i]=W[i]*D[Pre[i]]+C[i];
                    E[i]=W[i]*E[Pre[i]]+X[i];
                }
            }
            for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
                int i = *iter;
                if (flag[i] && !leaf[i]) T[i] = D[i]+E[i]/rho;
                int s1 = Suc1[i];
                int s2 = Suc2[i];
                if (!flag[s1]) T[s1]=T[i];
                if (!flag[s2]) T[s2]=T[i];
            }   
            //lambda = rho*Dlambda+Elambda
            double Dlambda=0;
            double Elambda=1./2.;
            if (flag[r] && flag[pr]){
                if (leaf[r]) Dlambda+=T[r]/(2*br);
                else {
                    Dlambda+=D[r]/(2*br);
                    Elambda+=E[r]/(2*br);
                }
                if (leaf[pr]) Dlambda-=T[pr]/(2*br);
                else {
                    Dlambda-=D[pr]/(2*br);
                    Elambda-=E[pr]/(2*br);
                }
            }
            else if (flag[r] && !flag[pr]){
                if (leaf[r]) Dlambda+=T[r]/(2*br);
                else {
                    Dlambda+=D[r]/(2*br);
                    Elambda+=E[r]/(2*br);
                }
                if (leaf[0]) Dlambda-=T[0]/(2*br);
                else {
                    Dlambda-=D[0]/(2*br);
                    Elambda-=E[0]/(2*br);
                }
            }
            else if (!flag[r] && flag[pr]){
                if (leaf[0]) Dlambda+=T[0]/(2*br);
                else {
                    Dlambda+=D[0]/(2*br);
                    Elambda+=E[0]/(2*br);
                }
                if (leaf[pr]) Dlambda-=T[pr]/(2*br);
                else {
                    Dlambda-=D[pr]/(2*br);
                    Elambda-=E[pr]/(2*br);
                }
            }
            lambda=Dlambda*rho+Elambda;
            B[r]=br*lambda;
            B[pr]=br-br*lambda;
            delete[] W;
            delete[] C;
            delete[] X;
            delete[] D;
            delete[] E;
            phi1 = phi(n,m,P,B,V,T,rho);
        }
        int *as= new int[m+3];
        for (int i=0;i<=m+2;i++) as[i]=-1;
        int count = 0;
        for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
            as[*iter]=count;
            count++;
        }
        bool *bl= new bool[m+1];
        for (int i=0;i<=m;i++) bl[i]=false;
        for (int y=0;y<ls.size();y++){
            while (!feuilles[y].empty()){
                int i=feuilles[y].top();
                if (!flag[i]){
                    feuilles[y].pop();
                    int p=P[i];
                    int s1=Suc1[i];
                    int s2=Suc2[i];
                    if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                            ldLagrange[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                        if (flag[s1] && flag[s2]) bl[i]=true;
                        if (!flag[s1] && bl[s1] && flag[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]];bl[i]=true;}
                        if (!flag[s2] && bl[s2] && flag[s1]) {ldLagrange[as[i]]+=ldLagrange[as[s2]];bl[i]=true;}
                        if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]]+ldLagrange[as[s2]];bl[i]=true;}
                    }
                    else{
                        int j=Suc1[p];
                        if (j==i)
                            j=Suc2[p];
                        ldLagrange[as[i]]=-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i]-2*rho*(B[j]-rho*T[j]+rho*T[p])/V[j];
                        if (P[p]!=-1)
                            ldLagrange[as[i]]+=2*rho*(B[p]-rho*T[p]+rho*T[P[p]])/V[p];
                        if (flag[j] && flag[p]) bl[i]=true;
                        if (flag[p] && !flag[j] && bl[j]) {ldLagrange[as[i]]-=ldLagrange[as[j]];bl[i]=true;}
                        if (!flag[p] && bl[p] && flag[j]) {ldLagrange[as[i]]+=ldLagrange[as[p]];bl[i]=true;}
                        if (!flag[j] && bl[j] && !flag[p] && bl[p]){ldLagrange[as[i]]+=ldLagrange[as[p]]-ldLagrange[as[j]];bl[i]=true;}
                    }
                }
            }
        }
        for (int y=0;y<top.size();y++){
            for (list<int>::iterator iter = internal[y].begin();iter!=internal[y].end();iter++){
                int i= *iter;
                int p = P[i];
                int s1=Suc1[i];
                int s2=Suc2[i];
                if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                        ldLagrange[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                    if (flag[s1] && flag[s2]) bl[i]=true;
                    if (!flag[s1] && bl[s1] && flag[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]];bl[i]=true;}
                    if (!flag[s2] && bl[s2] && flag[s1]) {ldLagrange[as[i]]+=ldLagrange[as[s2]];bl[i]=true;}
                    if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {ldLagrange[as[i]]+=ldLagrange[as[s1]]+ldLagrange[as[s2]];bl[i]=true;}
                }
            }
        }
        if (!flag[m+1]){
            ldLagrange[as[m+1]]=2*br*(br*lambda-rho*T[r]+rho*T[0])-2*br*(br-br*lambda-rho*T[pr]+rho*T[0]);
        }
        if (!flag[m+2]){
            ldLagrange[as[m+2]]=-2*br*(br*lambda-rho*T[r]+rho*T[0])+2*br*(br-br*lambda-rho*T[pr]+rho*T[0]);
        }
        for (int i=0;i<active_set.size();i++){
            if (myabs(ldLagrange[i])<(10e-15))
                ldLagrange[i]=0;
            ld.push_back(ldLagrange[i]);
        }
        delete[] SucL;
        delete[] SucI;
        delete[] feuilles;
        delete[] internal;
        delete[] Pre;
        delete[] leaf;
        delete[] as;
        delete[] bl;
        delete[] ldLagrange;
        return ld;
    }
}

double with_constraint_active_set_lambda_rate(double rho,int n,int m,double br,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &lambda){
    list<int> active_set;
    double phi1;
    phi1 = starting_point_lambda_rate(rho,n,m,br,P,Suc1,Suc2,B,V,T,lambda,active_set);
    double *T_new= new double[m+1];
    double *dir= new double[n];
    double alpha;
    for (int i=n;i<=m;i++) T_new[i]=T[i];
    bool *flag= new bool[m+3];
    double lambda_new;
    list<double> ldLagrange = without_constraint_active_set_lambda_rate(rho,n,m,br,P,flag,Suc1,Suc2,B,V,T_new,lambda_new,active_set,phi1);
    while (!conditions_lambda(m,ldLagrange,P,T_new,lambda_new)){
        for (int i=0;i<n;i++)  dir[i]=T_new[i]-T[i];
        double dirlambda = lambda_new-lambda;
        alpha=1;
        int as=0;
        double a;
        if (dirlambda<0){
            a = lambda/(-dirlambda);
            if (a<alpha){
                alpha=a;
                as = m+1;
            }
        }
        if (dirlambda>0){
            a = (1-lambda)/dirlambda;
            if (a<alpha){
                alpha=a;
                as=m+2;
            }
        }
        for (int i=n;i<=m;i++){
            if (flag[i]){
                int p = P[i];
                if (dir[p]>0){
                    a = (T[i]-T[p])/dir[p];
                    if (a<alpha){
                        alpha = a;
                        as = i;
                    }
                }
            }
        }
        for (int i=0;i<n;i++){
            if (flag[i]){
                int p = P[i];
                if (p!=-1){
                    if (dir[p]>dir[i]){
                        a = (T[i]-T[p])/(dir[p]-dir[i]);
                        if (a<alpha){
                            alpha = a;
                            as = i;
                        }
                    }
                }
            }
        }
        for (int i=0;i<n;i++) T[i]=T[i]+alpha*dir[i];
        lambda=lambda+alpha*dirlambda;
        int asrm= remove_ne_lambda(ldLagrange,active_set);
        if (asrm!=-1) active_set.remove(asrm);
        if (as!=0) active_set.push_back(as);
        ldLagrange = without_constraint_active_set_lambda_rate(rho,n,m,br,P,flag,Suc1,Suc2,B,V,T_new,lambda_new,active_set,phi1);
    }
    for (int i=0;i<n;i++) T[i]=T_new[i];
    lambda=lambda_new;
    delete[] T_new;
    delete[] dir;
    delete[] flag;
    return phi1;
}



double estimate_root_with_constraint_local_rate_rooted(double rho,int n,int m,int* &Suc1,int* &Suc2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda){
    //P: rooted tree, recherche la nouvelle racine autour de l'ancien racine.
    //estimate the root locally with constraint (QPD), the rate is known
    cout<<"Re-estimating the root with constraints around the given root ... ";
    double phi1=-1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    int s1=Suc1[0];
    int s2=Suc2[0];
    double br=reroot_rootedtree(m,s1,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    double ld;
    //cout<<"Optimizing the root position on the original branch "<<s1<<" ... ";
    cv[s1]=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
    //printf("%.10f\n",cv[s1]);
    cv[s2]=cv[s1];
    lambda=ld;
    r=s1;
    phi1=cv[s1];
    list<int> next;
    if (s1<n){ 
        next.push_back(Suc1[s1]);
        next.push_back(Suc2[s1]);
    }    
    if (s2<n){ 
        next.push_back(Suc1[s2]);
        next.push_back(Suc2[s2]);
    }
    while (!next.empty()){          
        int i = next.back();
        br=reroot_rootedtree(m,i,P,s1,s2,B,P_new,B_new);         
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<i<<" ... ";
        cv[i]=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
        //printf("%.10f\n",cv[i]);
        if (cv[i]<cv[P[i]]){
            if (i<n){
                next.push_back(Suc1[i]);
                next.push_back(Suc2[i]);
            }
            if (cv[i]<phi1){
                phi1=cv[i];r=i;lambda=ld;
            }
        }            
        next.remove(i);          
    }
    if (r==Suc1[0] || r==Suc2[0]) cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] cv;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_with_constraint_fast_rate_rooted(double rho,int n,int m,int s01,int s02,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda){
    //P: rooted tree, oublier la racine, recherche la nouvelle racine sur toutes les branches
    /////////////estimate root on all branches with QPD algorithm for rooted tree////////////////////////////////////////////////////////////
    cout<<"Re-estimating the root with constraints on all branches using fast method."<<endl;
    double ld;
    double phi1=estimate_root_without_constraint_rate_rooted(rho,n,m,s01,s02,P,B,V,T,r,ld);
    cout<<"Re-estimating the root with constraints around the given root ... ";    
    //cout<<"Optimizing the root position on the branch "<<r<<" ... ";
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    int* P_new = new int[m+1];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    int* P_ref = new int[m+1];
    int* tab = new int[m+1];
    double br=reroot_rootedtree(m,r,P,s01,s02,B,P_new,B_new,P_ref,tab);
    double* cv = new double[m+1];
    for (int i=0;i<=m;i++) cv[i]=0;
    computeSuc(P_new, Suc1_new, Suc2_new, m+1, n);
    cv[r]=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
    //printf("%.10f\n",cv[r]);
    list<int> next;
    phi1=cv[r];
    lambda=ld;
    int* Suc1_ref = new int[n];
    int* Suc2_ref = new int[n];
    computeSuc(P_ref,Suc1_ref,Suc2_ref,m+1,n);
    int s1=Suc1_ref[0];
    int s2=Suc2_ref[0];
    if (Suc1_ref[0]<n){
        next.push_back(Suc1_ref[s1]);
        next.push_back(Suc2_ref[s1]);
    }
    if (Suc2_ref[0]<n){
        next.push_back(Suc1_ref[s2]);
        next.push_back(Suc2_ref[s2]);
    }
    while (!next.empty()){
        int i = next.back();
        int e = tab[i];
        br=reroot_rootedtree(m,e,P,s01,s02,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<e<<" ... ";
        cv[e]=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
        //printf("%.10f\n",cv[e]);        
        if (cv[e]<cv[tab[P_ref[i]]]){
            if (i<n){
                next.push_back(Suc1_ref[i]);
                next.push_back(Suc2_ref[i]);
            }
            if (cv[e]<phi1){
                phi1=cv[e];r=e;lambda=ld;
            }
        }
        next.remove(i);
    }
    if (r==s01 || r==s02)  cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] cv;
    delete[] P_new;
    delete[] P_ref;
    delete[] tab;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] Suc1_ref;
    delete[] Suc2_ref;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_with_constraint_rate_rooted(double rho,int n,int m,int s1,int s2,int* & P,double* & B,double* &V,double* & T,int &r,double &lambda){
    //P: rooted tree
    cout<<"Re-estimating the root using constrained mode on all branches ... ";
    double phi1=-1;
    int y=1;
    int *P_new= new int[m+1];
    int *Suc1_new= new int[n];
    int *Suc2_new= new int[n];
    double *B_new= new double[m+1];
    double *T_new= new double[m+1];
    double ld;
    double br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
    //printf("%.10f\n",phi1);
    lambda=ld;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_rootedtree(m,y,P,s1,s2,B,P_new,B_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=with_constraint_active_set_lambda_rate(rho,n,m,br,P_new,Suc1_new,Suc2_new,B_new,V,T,ld);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
        }
        y++;
    }
    if (r==s1 || r==s2)  cout<<"The new root is on the original branch."<<endl;
    else cout<<"The new root is on the branch "<<r<<endl;
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}
