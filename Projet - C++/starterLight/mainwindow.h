#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


#define ARETE 0
#define TRIANGLE 1
#define INFINI 0
#define PAS_INFINI 1
#define MARK 0
#define UNMARK 1

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value; Color faceShadingColor;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;


enum DisplayMode {Normal, TemperatureMap, ColorShading};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // Fonction utilitaire
    MyMesh generer_3points();
    MyMesh generer_6points();
    MyMesh generer_3points_aligne();
    MyMesh generation_aleatoire(int n);
    void creer_boite_englobante(MyMesh *_mesh);
    void transfo_edge(MyMesh * _mesh, FaceHandle f);
    void recolor_edges(MyMesh * _mesh);
    void color_vertex(MyMesh * _mesh, int i, MyMesh::Color c);
    void creer_triangle_englobant(MyMesh *_mesh);
    bool est_ilegal_edge(MyMesh * _mesh, EdgeHandle e);
    void edge_flip_algo(MyMesh* _mesh);
    bool aucun_point_inf(MyMesh *_mesh, EdgeHandle e);
    int n_voisin_inf(MyMesh *_mesh, EdgeHandle e);
    std::vector<EdgeHandle> recup_voisin(MyMesh *_mesh, EdgeHandle e);
    // Geometrie
    float produit_scaliare(MyMesh::Point p, MyMesh::Point q);
    float produit_vectoriel(MyMesh::Point p, MyMesh::Point q);
    int determinant(MyMesh::Point A, MyMesh::Point B, MyMesh::Point C);
    bool est_dans_triangle(MyMesh * _mesh,VertexHandle p, FaceHandle fh);
    bool est_dans_arete(VertexHandle p, HalfedgeHandle heh);
    void separation_3_triangle(MyMesh *_mesh, VertexHandle p, FaceHandle f);
    void separation_4_triangle(MyMesh *_mesh, VertexHandle p, EdgeHandle e);
    void Hermeline(MyMesh * _mesh, int id_point);
    bool crit_boule_vide(MyMesh * _mesh, FaceHandle T, VertexHandle P);
    float rayon_cercle_circ(MyMesh::Point A, MyMesh::Point B, MyMesh::Point C);
    MyMesh::Point intersection(MyMesh::Point m1, MyMesh::Point v1, MyMesh::Point m2, MyMesh::Point v2);
    MyMesh::Point rotation(MyMesh::Point A, MyMesh::Point B);
    float distance(MyMesh::Point A, MyMesh::Point B);
    MyMesh::Point milieu(MyMesh::Point A, MyMesh::Point B);




    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:


    void on_horizontalSlider_valueChanged(int value);

    void on_horizontalSlider_2_valueChanged(int value);

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_spinBox_valueChanged(int arg1);

    void on_pushButton_6_clicked();

    void on_spinBox_2_valueChanged(int arg1);

    void on_pushButton_clicked();

private:
    bool modevoisinage;

    MyMesh mesh;
    MyMesh cloud;

    int nb_points;
    OpenMesh::VPropHandleT<int> Etat;
    OpenMesh::VPropHandleT<int> EdgeSupport;
    OpenMesh::VPropHandleT<int> TriangleSupport;
    OpenMesh::VPropHandleT<int> Inf;
    OpenMesh::FPropHandleT<bool> DEL;
    OpenMesh::EPropHandleT<bool> Mark;
    int id_point;
    bool separe_3;
    bool separe_4;
    int n_vertices;
    int id_edge;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
