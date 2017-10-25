 #include <iostream>
 #include <omp.h>
// #include <vector>
// #include <complex>
// #include <algorithm>
// #include <functional>
// #include <numeric>
// #include <ctime>

#include <Eigen/Eigen>
#include <omp.h>
#include <fstream>


using namespace Eigen;
 using namespace std;


void Meshrectangular(int dx ,int dy ,int nx, int ny , MatrixXd &Node, MatrixXd &Element);
void bcdof(VectorXd node, int dim, VectorXd &index);
void cal_BTDB (MatrixXd coord,double E, double nu,MatrixXd &B_T_D_B);
void maplocal(const VectorXd &index, const VectorXd &u_global , VectorXd & u_local);
void mapglobal(const VectorXd &index, VectorXd &u_global , VectorXd u_local);


int main() {
//  clock_t startTime = clock();
    // FILE *fp;
    // fp = fopen("example.dat","wb");

//cout<<mesh.active_local_elements_begin()<<"\n******\n";


// Domain Size
  double x_min = 0.0;
  double x_max = 100.0e3;
  double y_min = 0.0;
  double y_max = 50.0e3;
  int dim = 2.0;
  double dx = 50;
  double dy = 50;
   double nx = (x_max-x_min)/dx;
   double ny = (y_max-y_min)/dy;
  //double nx = 1;
  //double ny = 1;

  MatrixXd Node_top = MatrixXd::Zero((nx+1)*(ny+1),2);
  MatrixXd Element_top = MatrixXd::Zero(nx*ny,4);
  MatrixXd Node_bot = MatrixXd::Zero((nx+1)*(ny+1),2);
  MatrixXd Element_bot = MatrixXd::Zero(nx*ny,4);
  Meshrectangular(dx,dy,nx,ny,Node_top,Element_top);
  Meshrectangular(dx,dy,nx,ny,Node_bot,Element_bot);

  VectorXd x = VectorXd::LinSpaced(nx+1,x_min,x_max);
  VectorXd y = VectorXd::LinSpaced(ny+1,y_min,y_max);

  int station_neg_12_x = (((x_max-x_min)/2-12.0e3)/dx)+1;
  int station_neg_12_y = (0.0+3.0e3/dy)+1;
  int node_neg_12 = (station_neg_12_y-1)*(nx+1)+station_neg_12_x-1;
  int index_neg_12_x = node_neg_12*2;
  int index_neg_12_y = node_neg_12*2+1;

  int station_pos_12_x = (((x_max-x_min)/2+12.0e3)/dx)+1;
  int station_pos_12_y = (0.0+3.0e3/dy)+1;
  int node_pos_12 = (station_pos_12_y-1)*(nx+1)+station_pos_12_x-1;
  int index_pos_12_x = node_pos_12*2;
  int index_pos_12_y = node_pos_12*2+1;
  cout<<index_neg_12_x<<"\t"<<index_neg_12_y<<"\t"<<index_pos_12_x<<"\t"<<index_pos_12_y<<endl;
  int n_nodes = Node_top.rows();
  int nel = Element_top.rows();
  // Material Property
  double density = 2670.0;
  double v_s =3.464e3;
  double v_p = 6.0e3;
  double G= pow(v_s,2)*density;
  double Lambda = pow(v_p,2)*density-2.0*G;
  double E  = G*(3.0*Lambda+2.0*G)/(Lambda+G);
  double nu = Lambda/(2.0*(Lambda+G));

  // Time step
  double alpha = 0.4;
  double dt = alpha*dx/v_p;
  // Reyleigh Damping
  double beta =0.1;
  double q = beta*dt;
  // number of time steps

  double time_run = 12.0;

  int numt = time_run/dt;
  //int numt = 3;
  VectorXd time_total = dt*VectorXd::LinSpaced(numt,1,numt);
  // Number of dof per nodes
  int Ndofn = 2;
  // Number of node per elements
  int Nnel = 4;
  // Element mass
  double M=density*dx*dy*1.0;
  // External force vector
  VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);

  // solution vector
  // solution vector
  // + half plane
  VectorXd u_n_pos = VectorXd::Zero(n_nodes*dim,1);
  VectorXd v_n_pos = VectorXd::Zero(n_nodes*dim,1);
  VectorXd u_new_pos = VectorXd::Zero(n_nodes*dim,1);
  VectorXd v_new_pos = VectorXd::Zero(n_nodes*dim,1);
  VectorXd a_n_pos = VectorXd::Zero(n_nodes*dim,1);

  // - half plane
  VectorXd u_n_neg = VectorXd::Zero(n_nodes*dim,1);
  VectorXd v_n_neg = VectorXd::Zero(n_nodes*dim,1);
  VectorXd u_new_neg = VectorXd::Zero(n_nodes*dim,1);
  VectorXd v_new_neg = VectorXd::Zero(n_nodes*dim,1);
  VectorXd a_n_neg = VectorXd::Zero(n_nodes*dim,1);


  // slip and slip rate on the fault
  VectorXd delt_u_n_x = VectorXd::Zero(nx+1,1);
  VectorXd delt_v_n_x = VectorXd::Zero(nx+1,1);
  VectorXd delt_u_n_y = VectorXd::Zero(nx+1,1);
  VectorXd delt_v_n_y = VectorXd::Zero(nx+1,1);

  // Intial shear traction on the fault
  VectorXd Tx_0 = VectorXd::Zero(nx+1,1);
  for (int i=0 ; i<=nx; i++)
  {
    if ((x(i)<=(x_max+x_min)/2+1.5e3)&&(x(i)>=(x_max+x_min)/2-1.5e3))
    {
      Tx_0(i) = 81.6e6;
    }
    else if ((x(i)<=(x_max+x_min)/2+7.5e3+1.5e3)&&(x(i)>=(x_max+x_min)/2+7.5e3-1.5e3))
    {
      Tx_0(i) = 62.0e6;
    }
    else if ((x(i)<=(x_max+x_min)/2-7.5e3+1.5e3)&&(x(i)>=(x_max+x_min)/2-7.5e3-1.5e3))
    {
      Tx_0(i) = 78.0e6;
    }
    else
    {
         Tx_0(i) = 70.0e6;

    }
  }
  VectorXd sigmabar = 120.e6*VectorXd::Ones(nx+1,1);
  VectorXd Ty_0= -sigmabar;
  VectorXd T_cx = VectorXd::Zero(nx+1,1);
  VectorXd T_cy = VectorXd::Zero(nx+1,1);

  double Dc = 0.4;
  double mu_d=0.525;
  VectorXd mu_s = VectorXd::Zero(nx+1,1);
  for (int i=0; i<=nx;i++)
  {
    if ((x(i)<=(x_max+x_min)/2+15e3)&&(x(i)>=(x_max+x_min)/2-15e3))
    {
      mu_s(i) = 0.677;
    }
    else
    {
      mu_s(i) = 10000.0;
    }
  }
  VectorXd tau_s = VectorXd::Zero(nx+1,1);


  VectorXd M_el_vec = M/4*VectorXd::Ones(Nnel*Ndofn,1);


    VectorXd index_el_pos = VectorXd::Zero(Ndofn*Nnel,1);
    VectorXd index_el_neg = VectorXd::Zero(Ndofn*Nnel,1);

   MatrixXd index_store_pos = MatrixXd::Zero(Ndofn*Nnel,nel);
   MatrixXd index_store_neg = MatrixXd::Zero(Ndofn*Nnel,nel);

 //bcdof(Element.row(2),dim,index);
 for (int i=0;i<nel;i++)
 {
   bcdof(Element_top.row(i),dim,index_el_pos);
   bcdof(Element_bot.row(i),dim,index_el_neg);
   index_store_pos.col(i) = index_el_pos;
   index_store_neg.col(i) = index_el_neg;
 }
 printf("finisheindexing\n");
   VectorXd top_surf_index = VectorXd::Zero(2*(nx+1),1);
   VectorXd bot_surf_index = VectorXd::Zero(2*(nx+1),1);
   top_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,0,nx*2);
   top_surf_index.tail(nx+1)=VectorXd::LinSpaced(nx+1,1,nx*2+1);

   bot_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,2*ny*(nx+1),2*(ny+1)*(nx+1)-2);
   bot_surf_index.tail(nx+1) = VectorXd::LinSpaced(nx+1,2*ny*(nx+1)+1,2*(ny+1)*(nx+1)-1);

  // cout<<index_store_pos<<"\n"<<"*************"<<"\n"<<index_store_neg<<"\n";
    VectorXd M_global=VectorXd::Zero(n_nodes*Ndofn,1);
     for (int i=0 ; i<nel;i++)
      {
        index_el_pos = index_store_pos.col(i);
        mapglobal(index_el_pos,M_global,M_el_vec);
      }

//     // VectorXd coord_x = Node.block(0 , 0 , 3,0);
//      //VectorXd coord_y = Node.block(0 , 1 , 3,0);
// //  MatrixXd coord = Node.block(0,0,4,2);
  MatrixXd coord = MatrixXd::Zero(4,2);
   VectorXd Element_0= Element_top.row(0);
   coord.row(0) = Node_top.row(Element_0(0));
   coord.row(1) = Node_top.row(Element_0(1));
   coord.row(2) = Node_top.row(Element_0(2));
   coord.row(3) = Node_top.row(Element_0(3));
  MatrixXd B_T_D_B = MatrixXd::Zero(32,8);
   cal_BTDB (coord,E,nu,B_T_D_B);
   MatrixXd ke = MatrixXd::Zero(8,8);
   int k = 0;
   for (int i =0;i<4;i++)
   {
     ke = ke+ B_T_D_B.block(k,0,8,8);
     k=k+8;
   }
    // cout<<ke<<"\n";
//  //  Time Integration
// // cout<<index_store.col(0)<<index_store.col(1)<<endl;
//
printf("ready to start\n");
double start= omp_get_wtime();

std::ofstream slip_x("results/slip.txt");
std::ofstream rate_x("results/slip_rate.txt");
std::ofstream shear_x("results/shear.txt");
std::ofstream neg_12_u_x("results/neg_12_u_x.txt");
std::ofstream neg_12_u_y("results/neg_12_u_y.txt");
std::ofstream neg_12_v_x("results/neg_12_v_x.txt");
std::ofstream neg_12_v_y("results/neg_12_v_y.txt");

std::ofstream pos_12_u_x("results/pos_12_u_x.txt");
std::ofstream pos_12_u_y("results/pos_12_u_y.txt");
std::ofstream pos_12_v_x("results/pos_12_v_x.txt");
std::ofstream pos_12_v_y("results/pos_12_v_y.txt");


 for (int j=0;j<numt;j++)
  {
    VectorXd fe_global_pos= VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd fe_global_neg= VectorXd::Zero(n_nodes*Ndofn,1);
// //   cout<<index_store<<"\n"<<"*********"<<"\n";
//
double start = omp_get_wtime();
#pragma omp parallel
{
  #pragma omp for
  for (int i=0;i<nel;i++)
  {
    //printf("Element_num=%d",i);
    VectorXd u_n_local_pos=VectorXd::Zero(8,1) ;
    VectorXd v_n_local_pos=VectorXd::Zero(8,1) ;
    VectorXd u_n_local_neg=VectorXd::Zero(8,1) ;
    VectorXd v_n_local_neg=VectorXd::Zero(8,1) ;

    VectorXd index_pos = index_store_pos.col(i);
    VectorXd index_neg = index_store_neg.col(i);

    maplocal(index_pos,u_n_pos,u_n_local_pos);
    maplocal(index_pos,v_n_pos,v_n_local_pos);
    maplocal(index_neg,u_n_neg,u_n_local_neg);
    maplocal(index_neg,v_n_neg,v_n_local_neg);
  //  VectorXd fe_int_pos = VectorXd::Zero(8,1);
  //  VectorXd fe_int_neg = VectorXd::Zero(8,1);
 //   cout<<v_n<<"\n"<<"*********"<<"\n";

     VectorXd fe_int_pos = ke*(u_n_local_pos+q*v_n_local_pos);
     VectorXd fe_int_neg = ke*(u_n_local_neg+q*v_n_local_neg);
   // cout<<fe_int<<"\n"<<"*********"<<"\n";

    mapglobal(index_pos,fe_global_pos,fe_int_pos);
    mapglobal(index_neg,fe_global_neg,fe_int_neg);

  }
}
    VectorXd M_pos = VectorXd::Zero(nx+1,1);
    M_pos(0) = M/4.0;
    M_pos.segment(1,nx)=M/2*VectorXd::Ones(nx,1);
    M_pos(nx)= M/4.0;
    VectorXd M_neg=M_pos;

    VectorXd fe_int_fault_pos = VectorXd::Zero(nx+1,1);
    VectorXd fe_int_fault_neg = VectorXd::Zero(nx+1,1);

    maplocal(top_surf_index,fe_global_pos,fe_int_fault_pos);
    maplocal(bot_surf_index,fe_global_neg,fe_int_fault_neg);


    VectorXd fe_int_pos_x=-fe_int_fault_pos.head(nx+1);
     VectorXd fe_int_neg_x = -fe_int_fault_neg.head(nx+1);
    VectorXd fe_int_pos_y = fe_int_fault_pos.tail(nx+1);
     VectorXd fe_int_neg_y = -fe_int_fault_neg.tail(nx+1);
  //   cout<<top_surf_index<<"***********\n"<<bot_surf_index;
     double a = dx*1.0;
    VectorXd Tx = (1.0/dt*M_neg.cwiseProduct(M_pos).cwiseProduct(delt_v_n_x)+\
                  (M_neg.cwiseProduct(fe_int_pos_x)-M_pos.cwiseProduct(fe_int_neg_x))).cwiseQuotient(a*(M_neg+M_pos))+Tx_0;
    VectorXd Ty = (1.0/dt*M_neg.cwiseProduct(M_pos).cwiseProduct(delt_v_n_y+1.0/dt*(delt_u_n_y))+(M_neg.cwiseProduct(fe_int_pos_y)-M_pos.cwiseProduct(fe_int_neg_y))).cwiseQuotient(a*(M_neg+M_pos))+Ty_0;


    for (int k =0;k<nx+1;k++)
    {
      if (Ty(k)<=0.0)
      {
        T_cy(k) = Ty(k);
      }
      else
      {
        T_cy(k) = 0.0;
      }
    }
    for (int k =0;k<nx+1;k++)
    {
      if (delt_u_n_x(k)<Dc)
      {
        tau_s(k) = (mu_s(k)-(mu_s(k)-mu_d)*delt_u_n_x(k)/Dc)*(-T_cy(k));
      }
      else
      {
        tau_s(k) = mu_d*(-T_cy(k));
      }
    }
    for (int k =0;k<nx+1;k++)
    {
      if (Tx(k)<=tau_s(k))
      {
        T_cx(k) = Tx(k);
      }
      else
      {
        T_cx(k) = tau_s(k);
      }
    }
    VectorXd F_total_pos = F_ext_global-fe_global_pos;
     mapglobal(top_surf_index.head(nx+1),F_total_pos,-a*(T_cx-Tx_0));
     mapglobal(top_surf_index.tail(nx+1),F_total_pos,-a*(T_cy-Ty_0));

     VectorXd F_total_neg = F_ext_global-fe_global_neg;
     mapglobal(bot_surf_index.head(nx+1),F_total_neg,a*(T_cx-Tx_0));
     mapglobal(bot_surf_index.tail(nx+1),F_total_neg,a*(T_cy-Ty_0));
//
//   //  cout<<F_total(131);
//   // centerial difference solve
  a_n_pos = F_total_pos.cwiseQuotient(M_global);
  v_new_pos = v_n_pos+dt*a_n_pos;
  u_new_pos = u_n_pos+dt*v_new_pos;

  a_n_neg = F_total_neg.cwiseQuotient(M_global);
  v_new_neg = v_n_neg+dt*a_n_neg;
  u_new_neg = u_n_neg+dt*v_new_neg;
//
   v_n_pos = v_new_pos;
   u_n_pos = u_new_pos;

   v_n_neg = v_new_neg;
   u_n_neg = u_new_neg;//

   VectorXd u_n_fault_x_pos = VectorXd::Zero(nx+1,1);
   VectorXd v_n_fault_x_pos = VectorXd::Zero(nx+1,1);
   VectorXd u_n_fault_y_pos = VectorXd::Zero(nx+1,1);
   VectorXd v_n_fault_y_pos = VectorXd::Zero(nx+1,1);
   // neg
   VectorXd u_n_fault_x_neg = VectorXd::Zero(nx+1,1);
   VectorXd v_n_fault_x_neg = VectorXd::Zero(nx+1,1);
   VectorXd u_n_fault_y_neg = VectorXd::Zero(nx+1,1);
   VectorXd v_n_fault_y_neg = VectorXd::Zero(nx+1,1);
//
   maplocal(top_surf_index.head(nx+1),u_n_pos,u_n_fault_x_pos);
   maplocal(top_surf_index.head(nx+1),v_n_pos,v_n_fault_x_pos);
   maplocal(top_surf_index.tail(nx+1),u_n_pos,u_n_fault_y_pos);
   maplocal(top_surf_index.tail(nx+1),v_n_pos,v_n_fault_y_pos);
//
  maplocal(bot_surf_index.head(nx+1),u_n_neg,u_n_fault_x_neg);
  maplocal(bot_surf_index.head(nx+1),v_n_neg,v_n_fault_x_neg);
  maplocal(bot_surf_index.tail(nx+1),u_n_neg,u_n_fault_y_neg);
  maplocal(bot_surf_index.tail(nx+1),v_n_neg,v_n_fault_y_neg);
//
  delt_u_n_x = u_n_fault_x_pos-u_n_fault_x_neg;
  delt_u_n_y = u_n_fault_y_pos-u_n_fault_y_neg;
//
  delt_v_n_x = v_n_fault_x_pos-v_n_fault_x_neg;
  delt_v_n_y = v_n_fault_y_pos-v_n_fault_y_neg;
//   delt_v_n_x = 2.0*v_n_fault_x;
//
 //cout<<delt_u_n_x<<"\n"<<"**********"<<"\n";
 printf("Simulation time = %f\n",time_total(j));

 // fprintf(fp,"%f\n",delt_u_n_x);
 // fwrite(delt_u_n_x,fp);
 // fclose(fp);
 slip_x << delt_u_n_x<<"\n";
 rate_x << delt_v_n_x<<"\n";
 shear_x << T_cx<<"\n";
 neg_12_u_x << u_n_pos(index_neg_12_x)<<"\n";
 //cout<<u_n_pos(index_neg_12_x)<<"\n";
 neg_12_u_y << u_n_pos(index_neg_12_y)<<"\n";
  neg_12_v_x << v_n_pos(index_neg_12_x)<<"\n";
  neg_12_v_y << v_n_pos(index_neg_12_y)<<"\n";
 //
 pos_12_u_x << u_n_pos(index_pos_12_x)<<"\n";
 pos_12_u_y << u_n_pos(index_pos_12_y)<<"\n";
 pos_12_v_x << v_n_pos(index_pos_12_x)<<"\n";
 pos_12_v_y << v_n_pos(index_pos_12_y)<<"\n";


  }
  double end = omp_get_wtime();
  printf("CPU Time = %f\n",end-start);

//  //cout<<index_store<<"\n";
  return 0;

}
void Meshrectangular(int dx ,int dy ,int nx, int ny , MatrixXd &Node ,MatrixXd &Element)
{
  for (int i =0; i<=ny;i++)
  {
    for (int j= 0; j<=nx; j++)
    {
      Node((nx+1)*i+j,0)=j*dx;
      Node((nx+1)*i+j,1)=i*dy;
    }
  }
  for (int i =0; i<ny;i++)
  {
    for (int j= 0; j<nx; j++)
    {
      double Element_temp = nx*i+j;
      Element.row(nx*i+j) << Element_temp+i ,Element_temp+i+1 , Element_temp+i+nx+2, Element_temp+i+nx+1;
    }
  }
}
void bcdof(VectorXd node, int dim, VectorXd &index)
{
  // index.head(4)= dim*(node);
  // index.tail(4)= dim*(node)+VectorXd::Ones(4,1);
  index(0)=dim*(node(0));
  index(1)=dim*(node(0))+1;
  index(2)=dim*(node(1));
  index(3)=dim*(node(1))+1;
  index(4)=dim*(node(2));
  index(5)=dim*(node(2))+1;
  index(6)=dim*(node(3));
  index(7)=dim*(node(3))+1;
}
void cal_BTDB (MatrixXd coord ,double E, double nu, MatrixXd &B_T_D_B)
{
 MatrixXd D=MatrixXd::Zero(3,3);
  D << 1-nu, nu , 0,
      nu,   1-nu, 0,
      0,   0 ,  (1-2*nu)/2;
  D = E/((1+nu)*(1-2*nu))*D;
   VectorXd xi_gp=VectorXd::Zero(2,1);
   VectorXd eta_gp=VectorXd::Zero(2,1);
   VectorXd w_gp=VectorXd::Zero(2,1);
   xi_gp<<-1.0/sqrt(3),1.0/sqrt(3);
   eta_gp<<-1.0/sqrt(3),1.0/sqrt(3);
   w_gp<<1.0,1.0;
   int k = 0;
  for (int i=0; i<2;i++)
  {
    for (int j=0;j<2;j++)
    {
      double xi = xi_gp(i);
      double eta = eta_gp(j);
      double wgp = w_gp(j);
      // Jacobian Matrix (2x2)
      MatrixXd dNdxi=MatrixXd::Zero(2,4);
       dNdxi << eta-1.0, 1.0-eta, 1.0+eta, -eta-1.0,
               xi-1.0 , -xi-1.0, xi+1.0, 1.0-xi;
       dNdxi=dNdxi/4.0;
      //MatrixXd Jac=MatrixXd::Zero(2,2);
       MatrixXd Jac=MatrixXd::Zero(2,2);
       Jac = dNdxi*coord;
      //cout << "dNdxi:\n" << dNdxi << "\ncoord:\n" << coord << "\n\n";
      //cout<<Jac<<"\n\n";
   // Here is the problem
      //  Jac << 0.5,0,
      //         0 , 0.5;


      double detJ = Jac.determinant();
      MatrixXd dNdx = MatrixXd::Zero(2,4);
       dNdx = Jac.inverse()*dNdxi;
      MatrixXd B= MatrixXd::Zero(3,8);
      B<< dNdx(0,0),         0, dNdx(0,1),         0, dNdx(0,2),         0,  dNdx(0,3),         0,
                  0, dNdx(1,0),         0, dNdx(1,1),        0 , dNdx(1,2),          0,  dNdx(1,3),
          dNdx(1,0), dNdx(0,0), dNdx(1,1), dNdx(0,1), dNdx(1,2), dNdx(0,2),  dNdx(1,3),  dNdx(0,3);
      B_T_D_B.block(k,0,8,8)= B.transpose()*D*B*detJ*wgp;

      k=k+8;


    }
  }
}
void maplocal(const VectorXd &index, const VectorXd &u_global , VectorXd & u_local)
{
  for (int i=0;i<u_local.size();i++)
  {
    u_local(i)=u_global(index(i));
  }
}

void mapglobal(const VectorXd &index, VectorXd &u_global , VectorXd u_local)
{
  for (int i=0;i<u_local.size();i++)
  {
  //  u_local(i)=u_global(index(i));
    u_global(index(i))=u_global(index(i))+u_local(i);
  }
}
