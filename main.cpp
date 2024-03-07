#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <vector>
#include <sstream>
#include <time.h>
#include <armadillo>


#include<string>

using namespace std ;
using namespace arma ;

struct Courbe {
    mat points_de_controle; // Matrice des points de contrôle
    vec noeuds; // Vecteur des nœuds
    int degre; // degre de la B-spline

    Courbe(const mat& P, const vec& t, int p) : points_de_controle(P), noeuds(t), degre(p) {}
};

struct RepereFrenet {
    vec tangente;
    vec normale;
    vec binormale;
};

struct Point {
    double x;
    double y;
    double z;
};

struct Surface {
    mat points_de_controle; // Matrice des points de contrôle
    vec noeuds_u,
        noeuds_v; // Vecteur des nœuds
    int degre_u,
        degre_v; 

    Surface() {}

    Surface(const mat& P, const vec& t_u, const vec& t_v, int p_u, int p_v) 
            : points_de_controle(P), noeuds_u(t_u), noeuds_v(t_v), 
              degre_u(p_u), degre_v(p_v) {}
};

Surface surface;

#define M 4
#define N 4

std::vector<vec> pointsSurface;

double P[M][N][3];



void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);

void calculerPointsSurface(const Surface& surface, std::vector<vec>& pointsSurface);
double N_i_d(double u,int i,int d,const vec& t);

vec Surface_Nurbs(double u, double v, const Surface& surface);
// void afficherSurfaceNurbs(const Surface& surface);
void initSurface(void);
void afficherSurfaceNurbs() ;

RepereFrenet repereFrenet(const Courbe& courbe);
// vec derivee_1(const Courbe& courbe) ;
// vec derivee_2(const Courbe& courbe) ;

void afficheRepereFrenet(const Courbe& courbe, double u);

//void reshape(int,int);
float t=.5 ;

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
bool ctrlKeyPressed = false;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0.;

// constantes pour les materieux
  float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
    float mat_ambient[] = {0.7f, 0.7f, 0.7f, 1.0f};
    float mat_ambient_color[] = {0.8f, 0.8f, 0.2f, 1.0f};
    float mat_diffuse[] = {0.1f, 0.5f, 0.8f, 1.0f};
    float mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    float no_shininess = 0.0f;
    float low_shininess = 5.0f;
    float high_shininess = 100.0f;
    float mat_emission[] = {0.3f, 0.2f, 0.2f, 0.0f};

mat points(2, 4, arma::fill::zeros);

mat P_2(M*N, 3);
int pointCount = 0;

void initOpenGl() 
{ 

//lumiere 

	glClearColor( .5, .5, 0.5, 0.0 );
 
	glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat l_pos[] = { 3.,3.5,3.0,1.0 };
  glLightfv(GL_LIGHT0,GL_POSITION,l_pos);

  glLightfv(GL_LIGHT0,GL_DIFFUSE,l_pos);
 glLightfv(GL_LIGHT0,GL_SPECULAR,l_pos);
 glEnable(GL_COLOR_MATERIAL);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
// glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
         gluPerspective(45.0f,(GLfloat)200/(GLfloat)200,0.1f,10.0f);
	glMatrixMode(GL_MODELVIEW);
      gluLookAt(0.,0.,4., 0.,0.,0., 0.,1.,0.);

}

//------------------------------------------------------

void afficheBezier(mat P0, mat P1, mat P2, mat P3) 
{



  mat m = {{-1.0, 3.0, -3.0, 1.0},
            {3, -6.0, 3, 0.0},
            {-3.0, 3, 0, 1},
            {1, 0.0, 0, 0}};
  glColor3f(1.0, 1.0, 1.0); // Couleur blanche

  glBegin(GL_LINE_STRIP);
  for(double t = 0.0; t <= 1.0; t += 0.01) 
  {
    mat T = {t * t * t, t * t, t, 1};
    mat P = T * m * join_cols(P0, P1, P2, P3);
    glVertex3f(P[0], P[1], P[2]);
  }
  glEnd();
}

void afficheCatmullRom(mat P0, mat P1, mat P2, mat P3) 
{

  float s = 0.5;
  mat m = {{-1*s, 2-s, s-2, s},
            {2*s, s-3, 3-2*s, -1*s},
            {-1*s, 0, s, 0},
            {0, 1.0, 0, 0}};

  glColor3f(0.0, 1.0, 1.0); // Couleur blanche

  glBegin(GL_LINE_STRIP);
  for(double t = 0.0; t <= 1.0; t += 0.01) 
  {
    mat T = {t * t * t, t * t, t, 1};
    mat P = T * m * join_cols(P0, P1, P2, P3);
    glVertex3f(P[0], P[1], P[2]);
  }
  glEnd();
}

double N_i_d(double u,int i,int d,const vec& t)
{
  if(d==0)
  {
    if ((u >= t[i]) && (u < t[i+1]))
    {
      return 1.0;
    }
    return 0.0;
  }


    double Numérateur1 = u - t[i];
    double Denominateur1 = t[i + d] - t[i];
    double Numérateur2 = t[i + d + 1] - u;
    double Denominateur2 = t[i + d + 1] - t[i + 1];
      if (Denominateur1 == 0.0 || Denominateur2 == 0.0)
    {
        return 0.0;
    }


    double N1 = N_i_d(u, i, d-1, t);
    double N2 = N_i_d(u, i+1, d-1, t);

  return (Numérateur1/Denominateur1) * N1 + (Numérateur2/Denominateur2) * N2;
}

vec Courbe_B_Spline(double u, const Courbe& courbe) {
    vec C = zeros<vec>(courbe.points_de_controle.n_cols);

    for (int i = 0; i <= courbe.points_de_controle.n_cols; ++i) {
        double Ni_d = N_i_d(u, i, courbe.degre, courbe.noeuds);
        C += Ni_d * courbe.points_de_controle.row(i).t();
    }

    return C;
}

void afficheBSpline(const Courbe& courbe) {
    glColor3f(1.0, 0.0, 0.0);

    glBegin(GL_LINE_STRIP);

    for (double u = courbe.noeuds(courbe.degre); u <= courbe.noeuds(courbe.points_de_controle.n_rows); u += 0.01) {
        vec bspline = Courbe_B_Spline(u, courbe);
        glVertex3f(bspline(0), bspline(1), bspline(2));
    }

    glEnd();
}


// vec derivee_1(const Courbe& courbe) {
//     mat points = courbe.points_de_contrôle;
//     int n = points.n_cols;
//     vec derivee(n);

//     for (int i = 1; i < n - 1; ++i) {
//                 derivee(i) = (points(1, i + 1) - points(1, i - 1)) / (points(0, i + 1) - points(0, i - 1)); 
//     }

//     derivee(0) = (points(1, 1) - points(1, 0)) / (points(0, 1) - points(0, 0));
//     derivee(n - 1) = (points(1, n - 1) - points(1, n - 2)) / (points(0, n - 1) - points(0, n - 2));

//     return derivee;
// }

// vec derivee_2(const Courbe& courbe) {
//     mat points = courbe.points_de_contrôle;
//     int n = points.n_cols;
//     vec derivee(n);

//     for (int i = 1; i < n - 1; ++i) {
//         derivee(i) = (points(1, i + 1) - 2 * points(1, i) + points(1, i - 1)) 
//                     / ((points(0, i + 1) - points(0, i)) * (points(0, i) - points(0, i - 1))); 
//     }

//     derivee(0) = (points(1, 2) - 2 * points(1, 1) + points(1, 0)) 
//                 / ((points(0, 2) - points(0, 1)) * (points(0, 1) - points(0, 0)));
    
//     derivee(n - 1) = (points(1, n - 3) - 2 * points(1, n - 2) + points(1, n - 1)) 
//                     / ((points(0, n - 2) - points(0, n - 3)) * (points(0, n - 1) - points(0, n - 2)));

//     return derivee;
// }

vec calcule_derivee_c(const Courbe& courbe, double u)
{
    double h = 1e-6; //le pas
    vec C_plus_h = Courbe_B_Spline(u + h, courbe);
    vec C_minus_h = Courbe_B_Spline(u - h, courbe);

    return  (C_plus_h - C_minus_h) / (2 * h);
}

vec CalculerTangente( const Courbe& courbe, double u) {
    vec derivee_c=calcule_derivee_c(courbe, u);
    auto tangente = derivee_c / norm(derivee_c);

    return tangente;
}

vec calcule_derivee_c_2(const Courbe& courbe, double u)
 {
   double h = 1e-6;
    vec derivee_c = calcule_derivee_c(courbe, u);
    
    vec T = CalculerTangente(courbe, u);
    vec T_plus_h = CalculerTangente(courbe, u + h);
    vec T_minus_h = CalculerTangente(courbe, u - h);

    vec derivee_c_2 = (T_plus_h - 2 * T + T_minus_h) / (h * h);
    return derivee_c_2;
}

vec CalculerNormale(const Courbe& courbe, double u) {
    double h = 1e-6;  //le pas
    
    vec derivee_c =calcule_derivee_c(courbe, u);  
    vec derivee_c_2 = calcule_derivee_c_2(courbe, u);

    vec c_seconde_t = dot(derivee_c, derivee_c_2) * derivee_c / pow(norm(derivee_c), 2);
    vec c_seconde_n = derivee_c_2 - c_seconde_t;

   double norme_c_seconde_n = norm(cross(derivee_c_2, derivee_c)) / norm(derivee_c);
    vec normale = c_seconde_n / norme_c_seconde_n;

    if (dot(derivee_c_2, normale) < 0) {
        normale = -normale;
    }

    return normale;

    //     vec n =  derivee_2 / norm(derivee_2);
    //     return n;


}

double CalculerRayonCourbure(const Courbe& courbe, double u) {

    // Calcul du rayon de courbure R
    
    vec derivee_c = calcule_derivee_c(courbe, u);
    vec derivee_c_2 = calcule_derivee_c_2(courbe, u);

    double K = norm(cross(derivee_c, derivee_c_2)) / pow(norm(derivee_c), 3);
    double R = 1 / K;

    return R;
}


RepereFrenet repereFrenet(const Courbe& courbe, double u) {
    vec derivee_premiere = CalculerTangente(courbe, u);
    vec derivee_seconde = CalculerNormale(courbe, u);

    std::cout << "Dérivée première : " << derivee_premiere.t() << std::endl;
        std::cout << "Dérivée seconde : " << derivee_seconde.t() << std::endl;
    
    // tangente 
    vec tangente = normalise(derivee_premiere);

    //normal
    vec normale = normalise(derivee_seconde);

    //  binormal (produit vectoriel entre tangente et normal)
    vec binormale = normalise(cross(tangente, normale));

    return {tangente, normale, binormale};
}

void afficheRepereFrenet(const Courbe& courbe, double u) {
    RepereFrenet repere = repereFrenet(courbe,  u);

    vec position = Courbe_B_Spline(u, courbe);
    vec tangente = repere.tangente;
    vec normal = repere.normale;
    vec binormal = repere.binormale;

    // longueur des axes 
    double length = 0.1;
    
    // std::cout << "Vecteur tangente : " << repere.tangente.t() << std::endl;
    // std::cout << "Vecteur normal : " << repere.normale.t() << std::endl;
    // std::cout << "Vecteur binormal : " << repere.binormale.t() << std::endl;


    //  la tangente rouge
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(position(0), position(1), position(2));
    glVertex3f(position(0) + length * tangente(0), position(1) + length * tangente(1), position(2) + length * tangente(2));
    glEnd();

    // la normal vert
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(position(0), position(1), position(2));
    glVertex3f(position(0) + length * normal(0), position(1) + length * normal(1), position(2) + length * normal(2));
    glEnd();

    //  binormal bleu
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(position(0), position(1), position(2));
    glVertex3f(position(0) + length * binormal(0), position(1) + length * binormal(1), position(2) + length * binormal(2));
    glEnd();

}

void afficherCercleOsculateur(const Courbe& courbe, double u) {
    // Calcul du repère de Frenet au point u
    RepereFrenet repere = repereFrenet(courbe, u);

    // Récupération des vecteurs du repère de Frenet
    vec position = Courbe_B_Spline(u, courbe);
    vec tangente = CalculerTangente(courbe, u);
    vec normale = CalculerNormale(courbe, u);
    vec binormale = cross(tangente, normale);

    // Calcul du rayon de courbure au point u
    double rayon = CalculerRayonCourbure(courbe, u);
    vec centre = position + (rayon * normale);

    // Dessiner le cercle osculateur
    glBegin(GL_LINE_LOOP);
    for (double angle = 0; angle < 2 * M_PI; angle += 0.01) {
        double x = centre(0) + rayon * (cos(angle) * tangente(0) + sin(angle) * normale(0));
        double y = centre(1) + rayon * (cos(angle) * tangente(1) + sin(angle) * normale(1));
        double z = centre(2) + rayon * (cos(angle) * tangente(2) + sin(angle) * normale(2));
        glVertex3d(x, y, z);
    }
    glEnd();
}

void afficherSurfaceCourbe(const Courbe& courbe) {
    for (double t = 0.0; t <= 1.0; t += 0.005) {
        // Récupération des vecteurs de Frenet au point t
        RepereFrenet repere = repereFrenet(courbe, t);
        vec position = Courbe_B_Spline(t, courbe);

        vec tangente = repere.tangente;
        vec normal = repere.normale;
        vec binormale = repere.binormale;

        // Calcul du rayon de courbure au point t
        double rayon = 0.2;

        // Calcul du centre du cercle osculateur
        vec centre = position + (rayon * normal);

        // Dessiner le cercle osculateur
        glBegin(GL_LINE_LOOP);
        for (double angle = 0; angle < 2 * M_PI; angle += 0.01) {
               double x = position(0) + rayon * (cos(angle) * normal(0) + sin(angle) * binormale(0));
            double y = position(1) + rayon * (cos(angle) * normal(1) + sin(angle) * binormale(1));
            double z = position(2) + rayon * (cos(angle) * normal(2) + sin(angle) * binormale(2));
            glVertex3d(x, y, z);
        }
        glEnd();
    }
}

mat construire_matrice_P(const Courbe& courbe_1, const Courbe& courbe_2) {
    mat P = join_rows(courbe_1.points_de_controle, courbe_2.points_de_controle);
    return P;
}


void initGrille()
{


    // Remplissez le tableau avec des points de contrôle
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j)  {
            // Les coordonnées x et z sont régulièrement espacées de 0.25
            double x = i * 0.25;
            double z = j * 0.25;

            // Les coordonnées y peuvent être définies selon vos besoins
             double y = -0.2 + ((double)rand() / RAND_MAX) * 0.4;

            // Créez un point de contrôle
            vec point = {x, y, z};

            // Ajoutez le point au tableau
            P[i][j][0] = x;
            P[i][j][1] = y;
            P[i][j][2] = z;    
            // printf("Point [%d][%d]: x = %f, y = %f, z = %f\n", i, j, P[i][j][0], P[i][j][1], P[i][j][2]);   
          }
    }
   

}



void afficheSurfaceNubs(const Courbe& courbe_1, const Courbe& courbe_2, vec u, vec v) {
    int m = M - 1; // degre de la première courbe
    int n = N - 1; // degre de la deuxième courbe
    // int p = courbe_1.degre; // degre de la première courbe
    // int q = courbe_2.degre; // degre de la deuxième courbe


    glColor3f(1.0f, 0.0f, 1.0f);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            glPushMatrix();
            glTranslatef(P[i][j][0], P[i][j][1], P[i][j][2]);
            glutSolidSphere(0.05, 20, 20); // Ajustez le rayon et la résolution selon vos besoins
            glPopMatrix();
        }
    }
  
  glBegin(GL_TRIANGLES);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < 2; ++k) {
                for (int l = 0; l < 2; ++l) {
                    vec surface_point = {0, 0, 0};
                    // Calcul de la surface NURBS pour chaque sommet du triangle
                    for (int p = 0; p <= 1; ++p) {
                        for (int q = 0; q <= 1; ++q) {
                            double u_val = u[i + p];
                            double v_val = v[j + q];

                            double Ni = N_i_d(u_val, i + p, m, u); // Fonctions de base pour la courbe en u
                            double Nj = N_i_d(v_val, j + q, n, v); // Fonctions de base pour la courbe en v
                            
                            cout << "Ni = " << Ni << endl;
                            cout << "Nj = " << Nj << endl;

                            surface_point(0) += Ni * Nj * P[i + p][j+q][0];
                            surface_point(1) += Ni * Nj * P[i + p][j + q][1];
                            surface_point(2) += Ni * Nj * P[i + p][j + q][2];
                        }
                    }
                    glVertex3f(surface_point(0), surface_point(1), surface_point(2));
                    cout << "surface_point : " << surface_point << endl;
                }
            }
        }
    }
    glEnd();

}


mat initMatrice(mat& P) {
    
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j)  {
            double x = i * 4;
            double z = j * 4;

            double y = -2 + ((double)rand() / RAND_MAX) * 4;

            P.at(i * N + j, 0) = x;
            P.at(i * N + j, 1) = y;
            P.at(i * N + j, 2) = z;
        }
    }


   
    return P;
}


void initSurface(void)
{
  mat P_1(M * N, 3);
  P_1 = initMatrice(P_1);
  P_2 = P_1;
  cout << P_1 << endl;

  // mat P_1 = {
  //   {0.0, 0.0, 0.0},
  //   {1.0, 0.0, 0.0},
  //   {2.0, 0.0, 0.0},
  //   {3.0, 0.0, 0.0},

  //   {0.0, 1.0, 1.0},
  //   {1.0, 1.0, 1.0},
  //   {2.0, 1.0, 1.0},
  //   {3.0, 1.0, 1.0},

  //   {0.0, 2.0, 2.0},
  //   {1.0, 2.0, 2.0},
  //   {2.0, 2.0, 2.0},
  //   {3.0, 2.0, 2.0},

  //   {0.0, 3.0, 3.0},
  //   {1.0, 3.0, 3.0},
  //   {2.0, 3.0, 3.0},
  //   {3.0, 3.0, 3.0}
  // };

    int n_u = 3; 
    int d_u = 3; 
    int m_u = n_u + d_u + 1;

    int n_v = 3; 
    int d_v = 2; 
    int m_v = n_v + d_v + 1; 

    vec u = linspace<vec>(0, 1, m_u + 1);
    vec v = linspace<vec>(0, 1, m_v + 1);

    // Create the surface
    surface = Surface(P_1, u, v, d_u, d_v);
    cout << "surface done" << endl;

    // afficherSurfaceNurbs(surface);

  calculerPointsSurface(surface, pointsSurface);


}

vec Surface_Nurbs(double u, double v, const Surface& surface) {
    vec C = arma::zeros<vec>(surface.points_de_controle.n_rows);
//  cout << C << endl;
    for (int i = 0; i < surface.points_de_controle.n_rows; ++i) {

        double N_i_u = N_i_d(u, i, surface.degre_u, surface.noeuds_u);
        // cout << "N_i_u[" << i << "] = " << N_i_u << endl;
        for (int j = 0; j < surface.points_de_controle.n_cols; ++j) {
            double N_j_v = N_i_d(v, j, surface.degre_v, surface.noeuds_v);
            // cout << "N_j_v[" << j << "] = " << N_j_v << endl;
            // cout << "surface.points_de_contrôle(" << i << ", " << j << ") = " << surface.points_de_controle(i, j) << endl;
            C(i)+= N_i_u * N_j_v * surface.points_de_controle(i, j);
        }
    }
        cout << "Surface_Nurbs done" << endl;
            // cout << C << endl;

    return C;

}

void calculerPointsSurface(const Surface& surface, std::vector<vec>& pointsSurface) {
    int num_rows = surface.points_de_controle.n_rows;
    int num_cols = surface.points_de_controle.n_cols;
    double start_u = surface.noeuds_u(0);
    double end_u = surface.noeuds_u(surface.degre_u);
    double start_v = surface.noeuds_v(0);
    double end_v = surface.noeuds_v(surface.degre_v);
    double delta_u = (end_u - start_u) / (num_rows - 1);
    double delta_v = (end_v - start_v) / (num_cols - 1);

    for (int i = 0; i < num_rows - 1; ++i) {
        for (int j = 0; j < num_cols - 1; ++j) {
            double u1 = start_u + i * delta_u;
            double u2 = u1 + delta_u;
            double v1 = start_v + j * delta_v;
            double v2 = v1 + delta_v;

            vec surface_point1 = Surface_Nurbs(u1, v1, surface);
            vec surface_point2 = Surface_Nurbs(u2, v1, surface);
            vec surface_point3 = Surface_Nurbs(u1, v2, surface);
            vec surface_point4 = Surface_Nurbs(u2, v2, surface);

            pointsSurface.push_back(surface_point1);
            pointsSurface.push_back(surface_point2);
            pointsSurface.push_back(surface_point3);

            pointsSurface.push_back(surface_point2);
            pointsSurface.push_back(surface_point3);
            pointsSurface.push_back(surface_point4);
        }
    }
}

void afficherSurfaceNurbs() {
    glColor3f(0.0, 0.0, 1.0);

    glBegin(GL_TRIANGLES);

    // std::vector<vec> pointsSurface;
    // calculerPointsSurface(surface, pointsSurface);

    for (const auto& surface_point : pointsSurface) {
        glVertex3f(surface_point(0), surface_point(1), surface_point(2));
    }

    glEnd();

    
    
    // glColor3f(0.0, 0.0, 1.0);
    // glBegin(GL_LINES);
    // for (int i = 0; i < M; ++i) {
    //     for (int j = 0; j < N; ++j)  {
    //         if (i < M - 1) {
    //             glVertex3f(P_2(i * N + j, 0), P_2(i * N + j, 1), P_2(i * N + j, 2));
    //             glVertex3f(P_2((i + 1) * N + j, 0), P_2((i + 1) * N + j, 1), P_2((i + 1) * N + j, 2));
    //         }
    //         if (j < N - 1) {
    //             glVertex3f(P_2(i * N + j, 0), P_2(i * N + j, 1), P_2(i * N + j, 2));
    //             glVertex3f(P_2(i * N + j + 1, 0), P_2(i * N + j + 1, 1), P_2(i * N + j + 1, 2));
    //         }
    //     }
    // }
    // glEnd();


    cout << "afficherSurfaceNurbs done" << endl;
}





void displayCourbe(void)
{
  mat P0 = {1.0, 0.0, 0.0};
  mat P1 = {0.5, 0.5, -0.5};
  mat P2 = {-0.5, 0.5, 0.5};
  mat P3 = {0, 0.5, 0.5};

  mat P4 = {-1.0, 0.0, 0.0};
  mat P5 = {1, 1, 1};
  mat P6 = {-1, 0.5, 0.5};
  mat P7 = {0.0, -0.5, 0.5};

  mat P_1 = join_cols(P0, P1, P2, P3);
  mat P_2 = join_cols(P4, P5, P6, P7);

  int n1 = 3; // Nombre de points de contrôle - 1
  int d1 = 3; // degre de la B-spline
  int m1 = n1 + d1 + 1; // Nombre de nœuds

  int n2 = 3; // Nombre de points de contrôle - 1
  int d2 = 2; // degre de la B-spline
  int m2 = n2 + d2 + 1; // Nombre de nœuds

  vec u = linspace<vec>(0, 1, m1+1);
  vec v = linspace<vec>(0, 1, m2+1);

  Courbe courbe_1(P_1, u, d1);
    Courbe courbe_2(P_2, u, d2);

  afficheSurfaceNubs(courbe_1, courbe_2, u, v);


    // afficheBSpline(courbe_2);
//   afficheRepereFrenet(courbe_1, t);
//   afficheRepereFrenet(courbe_2, t);
// afficherCercleOsculateur(courbe_2, t) ;

  // afficheBSpline(courbe_1);

// afficherSurfaceCourbe(courbe_1);

  // GLfloat p0[] = {1.0, 0.0, 0.0};
  // GLfloat p1[] = {0.5, 0.5, -0.5};
  // GLfloat p2[] = {-0.5, 0.5, 0.5};
  // GLfloat p3[] = {0, 0.5, 0.5};

  // glPointSize(10.0);
  // glBegin(GL_POINTS);

  //   // Premier point
  //   glColor3f(1.0, 0.0, 0.0);  // Couleur rouge
  //   glVertex3fv(p0);

  //   // Deuxième point
  //   glColor3f(0.0, 1.0, 0.0);  // Couleur verte
  //   glVertex3fv(p1);

  //   // Troisième point
  //   glColor3f(0.0, 0.0, 1.0);  // Couleur bleue
  //   glVertex3fv(p2);

  //   // Quatrième point
  //   glColor3f(1.0, 1.0, 0.0);  // Couleur jaune
  //   glVertex3fv(p3);

  // glEnd();

  // afficheBezier(P0, P1, P2, P3);
  // afficheCatmullRom(P0, P1, P2, P3);


  // afficheBSpline(n1, d1, P_1, u);
  // afficheBSpline(n2, d2, P_2, v);

}


int main(int argc,char **argv)
{
  srand(time(NULL));
initGrille();

  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(600,600);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(1.0);
	
	//ifs = new Ifs();
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);
  //-------------------------------

  initSurface();

  //-------------------------------
    initOpenGl() ;
//-------------------------------

/* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}
//------------------------------------------------------
void affiche_repere(void)
{
  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(1.,0.);
  glEnd(); 

	 glBegin(GL_LINES);
  glColor3f(0.0,1.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(0.,1.);
  glEnd(); 
   glBegin(GL_LINES);
  glColor3f(0.0,0.0,1.0);
  glVertex3f(0.,0.,0.);
  glVertex3f(0.,0.,1.);
  glEnd(); 
}

//-----------------------------------------------------



//------------------------------------------------------
void affichage(void)
{
	glMatrixMode(GL_MODELVIEW);
  /* effacement de l'image avec la couleur de fond */
//	glClear(GL_COLOR_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//       glClearDepth(10.0f);                         // 0 is near, >0 is far

        glPushMatrix();
	glTranslatef(0,0,cameraDistance);
	glRotatef(cameraAngleX,1.,0.,0.)	;
	glRotatef(cameraAngleY,0.,1.,0.);
	// affiche_repere();
  // displayCourbe();
  afficherSurfaceNurbs();
        glPopMatrix();
  /* on force l'affichage du resultat */

          glFlush();
  glutSwapBuffers();

}

//------------------------------------------------------


//------------------------------------------------------
void clavier(unsigned char touche,int x,int y)
{

  switch (touche)
    {
    case '+': //
      t+=.01;
       if (t > 1 ) t=1;
      glutPostRedisplay();
      break;
    case '-': //* ajustement du t
       t-=.01;
        if (t < 0 ) t=0;
      glutPostRedisplay();
      break;
    case 'f': //* affichage en mode fil de fer 
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
      case 'p': //* affichage du carre plein 
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
  case 's' : //* Affichage en mode sommets seuls 
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;

    case 'q' : //*la touche 'q' permet de quitter le programme 
      exit(0);
    }
    
}
void mouse(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }

    else if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseMiddleDown = true;
        }
        else if(state == GLUT_UP)
            mouseMiddleDown = false;
    }
}

void mouseMotion(int x, int y)
{

    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        cameraDistance += (y - mouseY) * 0.2f;
        mouseY = y;
    }

    glutPostRedisplay();
}

    
    
