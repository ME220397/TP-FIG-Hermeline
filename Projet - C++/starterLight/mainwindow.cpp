#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "QRandomGenerator"

MyMesh MainWindow::generer_3points(){
    MyMesh nube;

    nube.add_property(Etat, "Etat");
    nube.add_property(EdgeSupport, "EdgeSupport");
    nube.add_property(TriangleSupport, "EdgeSupport");
    MyMesh::Point a(0, 0, 0);
    MyMesh::Point b(3, 0, 0);
    MyMesh::Point c(0, 3, 0);
    MyMesh::Point in(1,1, 0);
    MyMesh::Point mid(2,0,0);

    nube.add_vertex(a);
    nube.add_vertex(b);
    nube.add_vertex(c);
    nube.add_vertex(in);
    nube.add_vertex(mid);

    VertexHandle v1 = nube.vertex_handle(0);
    VertexHandle v2 = nube.vertex_handle(1);
    VertexHandle v3 = nube.vertex_handle(2);
    VertexHandle I = nube.vertex_handle(3);
    VertexHandle M = nube.vertex_handle(4);

    nube.data(v1).thickness = 10;
    nube.data(v2).thickness = 10;
    nube.data(v3).thickness = 10;
    nube.data(I).thickness = 10;
    nube.data(M).thickness = 10;

    nube.set_color(v1, MyMesh::Color(0, 0, 255));
    nube.set_color(v2, MyMesh::Color(0, 0, 255));
    nube.set_color(v3, MyMesh::Color(0, 0, 255));

    nube.add_face(v1, v2, v3);

    EdgeHandle e1 = nube.edge_handle(0);
    EdgeHandle e2 = nube.edge_handle(1);
    EdgeHandle e3 = nube.edge_handle(2);

    nube.data(e1).thickness = 2;
    nube.data(e2).thickness = 2;
    nube.data(e3).thickness = 2;

    nube.set_color(e1, MyMesh::Color(0, 0, 0));
    nube.set_color(e2, MyMesh::Color(0, 0, 0));
    nube.set_color(e3, MyMesh::Color(0, 0, 0));

    assert(est_dans_triangle(&nube, I, nube.face_handle(0)) == true);
    assert(est_dans_triangle(&nube, M, nube.face_handle(0)) == true);


    return nube;
}

MyMesh MainWindow::generer_6points(){
    MyMesh nube;

    nube.add_property(Etat, "Etat");
    nube.add_property(EdgeSupport, "EdgeSupport");
    nube.add_property(TriangleSupport, "EdgeSupport");
    MyMesh::Point a(0, 0, 0);
    MyMesh::Point b(4, 0, 0);
    MyMesh::Point c(0, 4, 0);
    MyMesh::Point d(4, 4, 0);
    MyMesh::Point in(1,1, 0);
    MyMesh::Point mid(2,2,0);

    nube.add_vertex(a);
    nube.add_vertex(b);
    nube.add_vertex(c);
    nube.add_vertex(d);
    nube.add_vertex(in);
    nube.add_vertex(mid);

    VertexHandle v1 = nube.vertex_handle(0);
    VertexHandle v2 = nube.vertex_handle(1);
    VertexHandle v3 = nube.vertex_handle(2);
    VertexHandle v4 = nube.vertex_handle(3);
    VertexHandle I = nube.vertex_handle(4);
    VertexHandle M = nube.vertex_handle(5);

    nube.data(v1).thickness = 10;
    nube.data(v2).thickness = 10;
    nube.data(v3).thickness = 10;
    nube.data(v4).thickness = 10;
    nube.data(I).thickness = 10;
    nube.data(M).thickness = 10;

    nube.set_color(v1, MyMesh::Color(0, 0, 255));
    nube.set_color(v2, MyMesh::Color(0, 0, 255));
    nube.set_color(v3, MyMesh::Color(0, 0, 255));
    nube.set_color(v4, MyMesh::Color(0, 0, 255));

    nube.add_face(v1, v2, v3);
    nube.add_face(v4, v3, v2);

    for(MyMesh::EdgeIter e_i = nube.edges_begin(); e_i != nube.edges_end(); e_i++){
        nube.data(*e_i).thickness = 2;
        nube.set_color(*e_i, MyMesh::Color(0, 0, 0));
    }

    FaceHandle t1 = nube.face_handle(0);
    FaceHandle t2 = nube.face_handle(1);

    assert(est_dans_triangle(&nube, I, t1) == true);
    assert(est_dans_triangle(&nube, M, t1) == true);

    assert(est_dans_triangle(&nube, I, t2) == false);
    assert(est_dans_triangle(&nube, M, t2) == true);
    if(est_dans_triangle(&nube, I, t1)){
        if(nube.property(Etat, I) == TRIANGLE)
            if(separe_3)
                separation_3_triangle(&nube, I, t1);
    }

    if(est_dans_triangle(&nube, M, t1)){
        if(nube.property(Etat, M) == ARETE){
            if(separe_4){
                int id = nube.property(EdgeSupport, M);
                separation_4_triangle(&nube, M, nube.edge_handle(id));
            }
        }
    }
    recolor_edges(&nube);
    return nube;
}

MyMesh MainWindow::generation_aleatoire(int n){
    MyMesh nube;
    nube.add_property(Etat, "Etat");
    nube.add_property(EdgeSupport, "EdgeSupport");
    nube.add_property(TriangleSupport, "EdgeSupport");
    nube.add_property(Inf, "Inf");
    nube.add_property(DEL, "DEL");
    int min = -25;
    int max = 25;
    for (int i=0; i<n; i++) {
        int x = QRandomGenerator::global()->bounded(min, max);
        int y = QRandomGenerator::global()->bounded(min, max);

        nube.add_vertex(MyMesh::Point(x, y, 0));
    }

    for (MyMesh::VertexIter v = nube.vertices_begin(); v != nube.vertices_end(); v++) {
        nube.data(*v).thickness = 10;
        nube.set_color(*v, MyMesh::Color(0, 0, 255));
    }

    //creer_boite_englobante(&nube);
    creer_triangle_englobant(&nube);
    n_vertices = nube.n_vertices();
    id_point = 0;
    color_vertex(&nube, id_point, MyMesh::Color(255, 165, 0));
    return nube;
}

void colorer_point_boite(MyMesh *_mesh, VertexHandle v){
    _mesh->data(v).thickness = 10;
    _mesh->set_color(v, MyMesh::Color(255, 0, 0));
}

void MainWindow::creer_boite_englobante(MyMesh *_mesh){
    int minx = 40, miny = 40;
    int maxx = -40, maxy = -40;

    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); v_it++){
        MyMesh::Point p = _mesh->point(*v_it);
        if(minx > p[0])
          minx = p[0];
        if(miny > p[1])
          miny = p[1];
        if(maxx < p[0])
          maxx = p[0];
        if(maxy < p[1])
          maxy = p[1];
    }

    MyMesh::Point A(minx - 2, miny - 2, 0);
    MyMesh::Point B(maxx + 2, maxy + 2, 0);

    MyMesh::Point C(maxx + 2, miny - 2, 0);
    MyMesh::Point D(minx - 2, maxy + 2, 0);

    VertexHandle v1 = _mesh->add_vertex(A);
    VertexHandle v2 = _mesh->add_vertex(B);
    VertexHandle v3 = _mesh->add_vertex(C);
    VertexHandle v4 = _mesh->add_vertex(D);

    // Rajouter la propriété infini aux points
    _mesh->property(Inf, v1) = INFINI;
    _mesh->property(Inf, v2) = INFINI;
    _mesh->property(Inf, v3) = INFINI;
    _mesh->property(Inf, v4) = INFINI;

    colorer_point_boite(_mesh, v1);
    colorer_point_boite(_mesh, v2);
    colorer_point_boite(_mesh, v3);
    colorer_point_boite(_mesh, v4);

    FaceHandle t1 = _mesh->add_face(v1, v3, v2);
    FaceHandle t2 = _mesh->add_face(v1, v2, v4);

    transfo_edge(_mesh, t1);
    transfo_edge(_mesh, t2);
}

void MainWindow::creer_triangle_englobant(MyMesh *_mesh){
    int minx = 40, miny = 40;
    int maxx = -40, maxy = -40;

    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); v_it++){
        MyMesh::Point p = _mesh->point(*v_it);
        if(minx > p[0])
          minx = p[0];
        if(miny > p[1])
          miny = p[1];
        if(maxx < p[0])
          maxx = p[0];
        if(maxy < p[1])
          maxy = p[1];
    }

    MyMesh::Point A(minx - 60, miny - 20, 0);
    MyMesh::Point B(maxx + 60, miny - 20, 0);
    MyMesh::Point C((A[0] + B[0]) / 2, maxy + 60, 0);

    VertexHandle v1 = _mesh->add_vertex(A);
    VertexHandle v2 = _mesh->add_vertex(B);
    VertexHandle v3 = _mesh->add_vertex(C);

    // Rajouter la propriété infini aux points
    _mesh->property(Inf, v1) = INFINI;
    _mesh->property(Inf, v2) = INFINI;
    _mesh->property(Inf, v3) = INFINI;

    colorer_point_boite(_mesh, v1);
    colorer_point_boite(_mesh, v2);
    colorer_point_boite(_mesh, v3);

    FaceHandle t1 = _mesh->add_face(v1, v2, v3);

    transfo_edge(_mesh, t1);
}


// Geometrie

float MainWindow::produit_scaliare(MyMesh::Point p, MyMesh::Point q){
   float px, py, pz;
   float qx, qy, qz;
   px = p[0]; py = p[1]; pz = p[2];
   qx = q[0]; qy = q[1]; qz = q[2];

   return px*qx + py*qy + pz*qz;
}

void MainWindow::Hermeline(MyMesh *_mesh, int id_point){

    VertexHandle v = _mesh->vertex_handle(id_point);
    // On itere sur les points
    if(_mesh->property(Inf, v) != INFINI){
        for(MyMesh::FaceIter f_it = _mesh->faces_begin(); f_it!= _mesh->faces_end(); f_it++){
            FaceHandle fh = *f_it;
            if(fh.is_valid() && est_dans_triangle(_mesh, v, fh)){
                if(_mesh->property(Etat, v) == TRIANGLE){
                    separation_3_triangle(_mesh, v, fh);
                    _mesh->property(Etat, v) = -1;
                }
                if(_mesh->property(Etat, v) == ARETE){
                    int id = _mesh->property(EdgeSupport, v);
                    EdgeHandle e = _mesh->edge_handle(id);
                    separation_4_triangle(_mesh, v, e);
                    _mesh->property(Etat, v) = -1;
                }
                break;
            }
        }
    }

}

int MainWindow::determinant(MyMesh::Point A, MyMesh::Point B, MyMesh::Point P){
    // Un fois que l'on a recuperé les points du triangle
    // On verifie qu'un point P se situe a gauche de AB
    /*
                | Ax Bx Px |        det > 0 -> P à gauche
          det = | Ay By Py |        det < 0 -> P à droite
                | 1  1  1  |        det = 0 -> P inclu dans le segment

                  | By Py |      | Ax Py |      | Ay Px |
          det = Ax| 1  1  | - Bx | 1  1  | + Px | 1  1 |
     */
    int x = 0, y = 1;
    int det = A[x] * (B[y] - P[y]) - B[x] * (A[y] - P[y]) + P[x] * (A[y] - B[y]);

    return det;
}

void MainWindow::separation_3_triangle(MyMesh *_mesh,VertexHandle P, FaceHandle f){
    int old_n = _mesh->n_faces();
    // On recupére les point v0, v1 et v2 du triangle
    MyMesh::VertexHandle vh_A, vh_B, vh_C;
    MyMesh::Point A, B, C;
    MyMesh::FaceVertexIter v_it = _mesh->fv_iter(f);

    if(v_it.is_valid()){
        vh_A = *v_it;
        A = _mesh->point(vh_A);
    }
    v_it++;
    if(v_it.is_valid()){
        vh_B = *v_it;
        B = _mesh->point(vh_B);
    }
    v_it++;
    if(v_it.is_valid()){
        vh_C = *v_it;
        C = _mesh->point(vh_C);
    }

    // On supprime le triangle une fois les point recupéré
    _mesh->request_face_status();
    _mesh->request_edge_status();
    _mesh->request_vertex_status();
    _mesh->delete_face(f, false);
    _mesh->garbage_collection();

    // A -> v0, B -> v1, C -> v2
    int det_ABC = determinant(A, B, C);
    if(det_ABC > 0 ){
        _mesh->add_face(P, vh_A, vh_B);
        _mesh->add_face(P, vh_B, vh_C);
        _mesh->add_face(P, vh_C, vh_A);
    }
    else if(det_ABC < 0){
        _mesh->add_face(P, vh_B, vh_A);
        _mesh->add_face(P, vh_A, vh_C);
        _mesh->add_face(P, vh_C, vh_B);
    }
    // On s'assure que l'on à bien 2 triangle en plus
    int n = _mesh->n_faces();
    assert(n == old_n+2);
}
void MainWindow::transfo_edge(MyMesh * _mesh, FaceHandle f){
    for(MyMesh::FaceEdgeIter fe_it = _mesh->fe_iter(f); fe_it.is_valid(); fe_it++){
        EdgeHandle eh = *fe_it;
        _mesh->data(eh).thickness = 2;
        _mesh->set_color(eh, MyMesh::Color(0, 0, 0));
    }
}
void MainWindow::separation_4_triangle(MyMesh *_mesh ,VertexHandle P, EdgeHandle e){
    int n_old = _mesh->n_faces();
    //On recupere les HalfedgeHandle permettant de facilité la generation de nouveau triangle
    HalfedgeHandle he_t1 = _mesh->halfedge_handle(e, 0);
    HalfedgeHandle he_t2 = _mesh->halfedge_handle(e, 1);
    HalfedgeHandle he_t1_next = _mesh->next_halfedge_handle(he_t1);
    HalfedgeHandle he_t2_next = _mesh->next_halfedge_handle(he_t2);
    // On recupere les deux triangles que l'on va supprimer
    FaceHandle t1 = _mesh->face_handle(he_t1);
    FaceHandle t2 = _mesh->face_handle(he_t2);

    // Sommets des deux triangles
    VertexHandle vh_A, vh_B, vh_C, vh_D;
    MyMesh::Point A, B, C, D;

    // On commence par le triangle T1
    vh_A = _mesh->from_vertex_handle(he_t1);
    vh_B = _mesh->to_vertex_handle(he_t1);
    vh_C = _mesh->to_vertex_handle(he_t1_next);
    vh_D = _mesh->to_vertex_handle(he_t2_next);

    A = _mesh->point(vh_A);
    B = _mesh->point(vh_B);
    C = _mesh->point(vh_C);
    D = _mesh->point(vh_D);

    _mesh->request_face_status();
    _mesh->request_edge_status();
    _mesh->request_vertex_status();
    _mesh->delete_face(t1, false);
    _mesh->delete_face(t2, false);
    _mesh->garbage_collection();

    int det_ABC = determinant(A, B ,C);
    if(det_ABC > 0){
        // Creation de deux triangle a partir de ABC
        _mesh->add_face(P, vh_B , vh_C);
        _mesh->add_face(P, vh_C , vh_A);

        // Creation de deux triangle a partir de ABD
        _mesh->add_face(P, vh_A , vh_D);
        _mesh->add_face(P, vh_D , vh_B);
    }
    else if(det_ABC < 0){
        // Creation de deux triangle a partir de ABC
        _mesh->add_face(P, vh_B , vh_D);
        _mesh->add_face(P, vh_D , vh_A);

        // Creation de deux triangle a partir de ABD
        _mesh->add_face(P, vh_A , vh_C);
        _mesh->add_face(P, vh_C , vh_B);
    }

    // On supprime les deux triangles d'origine
    // On s'assure que l'on a bien 2 faces de plus qu'avant
    int n = _mesh->n_faces();
    qDebug() << "n = " << n << " n_old = " << n_old;
    assert(n == n_old+2);
}
bool MainWindow::est_dans_triangle(MyMesh * _mesh, VertexHandle M, FaceHandle fh){
    MyMesh::Point A, B, C;
    MyMesh::Point P = _mesh->point(M);
    HalfedgeHandle AB = _mesh->halfedge_handle(fh);
    VertexHandle v1 = _mesh->from_vertex_handle(AB);
    VertexHandle v2 = _mesh->to_vertex_handle(AB);
    HalfedgeHandle BC = _mesh->next_halfedge_handle(AB);
    VertexHandle v3 = _mesh->to_vertex_handle(BC);
    HalfedgeHandle CA = _mesh->next_halfedge_handle(BC);

    A = _mesh->point(v1);
    B = _mesh->point(v2);
    C = _mesh->point(v3);

    int x = 0, y = 1;


    // Segment A B

    float det1 = determinant(A, B, P);
    if(det1 < 0)
        return false;
    if(det1 == 0){
        _mesh->property(Etat, M) = ARETE;
        EdgeHandle e =  _mesh->edge_handle(AB);
        _mesh->property(EdgeSupport, M) = e.idx();
        return true;
    }

    // Segment B C
    float det2 = determinant(B, C, P);
    if(det2 < 0)
        return false;
    if(det2 == 0){
        _mesh->property(Etat, M) = ARETE;
        EdgeHandle e =  _mesh->edge_handle(BC);
        _mesh->property(EdgeSupport, M) = e.idx();
        return true;
    }

    // Segment C A
    float det3 = determinant(C, A, P);
    if(det3 < 0)
        return false;
    if(det3 == 0){
        _mesh->property(Etat, M) = ARETE;
        EdgeHandle e =  _mesh->edge_handle(CA);
        _mesh->property(EdgeSupport, M) = e.idx();
        return true;
    }

    // Si on arrive a ce points alors tout les det sont supérieur à 0
    _mesh->property(Etat, M) = TRIANGLE;
    _mesh->property(TriangleSupport, M) = fh.idx();
    return true;

}

MyMesh::Point milieu(MyMesh::Point A, MyMesh::Point B)
{
    MyMesh::Point I;
    I[0] = (A[0] + B[0])/2;
    I[1] = (A[1] + B[1])/2;
    return I;
}


float distance(MyMesh::Point A, MyMesh::Point B)
{
    float dist = sqrt(pow(A[0] - B[0], 2) + pow(A[1] - B[1], 2));
    return dist;
}


MyMesh::Point rotation(MyMesh::Point A, MyMesh::Point B)
{
    MyMesh::Point AB = B-A;
    MyMesh::Point mediatrice;
    mediatrice[0] = AB[0]*cos(90) - AB[1]*sin(90);
    mediatrice[1] = AB[0]*sin(90) + AB[1]*cos(90);
    return mediatrice;
}


MyMesh::Point intersection(MyMesh::Point m1, MyMesh::Point v1, MyMesh::Point m2, MyMesh::Point v2)
{
    MyMesh::Point inter;

    float t = (v2[0]*(m2[1]-m1[2]) + v2[1]*(m2[0]-m2[0]))/(v2[0]*v1[1] - v1[0]*v2[1]);
    inter[0] = m1[0] + t*v1[0];
    inter[1] = m1[1] + t*v1[1];
    return inter;
}


bool crit_boule_vide(MyMesh * _mesh, EdgeHandle eh)
{
    HalfedgeHandle heh = _mesh->halfedge_handle(eh,0);
    HalfedgeHandle heh_opp = _mesh->halfedge_handle(eh,1);
    //Initialisation
    VertexHandle Ah = _mesh->from_vertex_handle(heh);
    VertexHandle Bh = _mesh->to_vertex_handle(heh);
    MyMesh::Point A = _mesh->point(Ah);
    MyMesh::Point B = _mesh->point(Bh);
    //recuperation du milieu de [AB]
    MyMesh::Point milieu1 = milieu(A,B);
    //Recuperation de la mediatrice du segment AB
    MyMesh::Point mediatrice1 = rotation(A,B);

    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle Ch = _mesh->to_vertex_handle(heh);
    MyMesh::Point C = _mesh->point(Ch);
    //de meme pour le segment BC
    MyMesh::Point milieu2 = milieu(B,C);
    MyMesh::Point mediatrice2 = rotation(B,C);

    //recuperation du point du triangle oppose
    heh_opp = _mesh->next_halfedge_handle(heh_opp);
    VertexHandle Dh = _mesh->to_vertex_handle(heh_opp);
    MyMesh::Point D = _mesh->point(Dh);
    /*
    MyMesh::Point milieu3 = milieu(A,D);
    MyMesh::Point mediatrice3 = rotation(A,D);*/

    MyMesh::Point centre_circ = intersection(milieu1, mediatrice1, milieu2, mediatrice2);
    float rayon = distance(A, centre_circ);
    if(rayon >= distance(D, centre_circ))
        return true;

/*    centre_circ = intersection(milieu1, mediatrice1, milieu3, mediatrice3);
    rayon = distance(D, centre_circ);
    if(rayon >= distance(A, centre_circ))
        return true;*/

    return false;
}













// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, DisplayMode mode)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        qSort(values);

        float range = values.at(values.size()*0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::Normal)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    //ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    separe_3 = false;
    separe_4 = false;
    nb_points = 15;
    id_point = 0;
    n_vertices = 0;
    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_clicked()
{
    cloud = generer_3points();
    displayMesh(&cloud);
}

void MainWindow::on_horizontalSlider_valueChanged(int value)
{
    MyMesh * _mesh = &cloud;

    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        _mesh->data(*v).thickness = value;
    }
    displayMesh(&cloud);
}

void MainWindow::on_horizontalSlider_2_valueChanged(int value)
{
    MyMesh * _mesh = &cloud;

    for(MyMesh::EdgeIter e = _mesh->edges_begin(); e != _mesh->edges_end(); e++){
        _mesh->data(*e).thickness = value;
    }
    displayMesh(&cloud);
}

void MainWindow::on_pushButton_2_clicked()
{
    cloud = generer_6points();
    displayMesh(&cloud);
}

void MainWindow::on_pushButton_3_clicked()
{
    separe_3 = true;
    cloud = generer_6points();
    separe_3 = false;
    displayMesh(&cloud);
}

void MainWindow::on_pushButton_4_clicked()
{
    separe_4 = true;
    cloud = generer_6points();
    separe_4 = false;
    displayMesh(&cloud);
}

void MainWindow::on_pushButton_5_clicked()
{
    cloud = generation_aleatoire(nb_points);
    recolor_edges(&cloud);
    displayMesh(&cloud);
}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    nb_points = arg1;
}

void MainWindow::recolor_edges(MyMesh * _mesh){
    for(MyMesh::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); e_it++){
        _mesh->data(*e_it).thickness = 2;
        _mesh->set_color(*e_it, MyMesh::Color(0,0,0));
    }
}

void MainWindow::color_vertex(MyMesh * _mesh, int i, MyMesh::Color c){
    VertexHandle v = _mesh->vertex_handle(i);
    _mesh->set_color(v, c);
}
void MainWindow::on_pushButton_6_clicked()
{
    if(id_point < nb_points){
        Hermeline(&cloud, id_point);
        color_vertex(&cloud, id_point++, MyMesh::Color(0, 0, 255));
        recolor_edges(&cloud);
        color_vertex(&cloud, id_point, MyMesh::Color(255, 165, 0));
    }
    displayMesh(&cloud);
}
