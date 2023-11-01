#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;


// input processor : change typeSpec to nodes/r for 2 triangle case
void inputProcessorFunc(const std::string& inputFileName, vector<Vector3d>* Nodes, vector<vector<int>>* facenodes) {
    ifstream inputFile(inputFileName);
    string line, typeSpec;

    Nodes->clear();
    facenodes->clear();

    while (getline(inputFile, line)) {
        // Trim leading and trailing spaces
        line.erase(line.find_last_not_of(" \t") + 1);
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty()) {
            continue;
        } else if (line[0] == '*') {
            // Set typeSpec without leading spaces and convert to lowercase for case-insensitive comparison
            typeSpec = line.substr(1);
            transform(typeSpec.begin(), typeSpec.end(), typeSpec.begin(), ::tolower);
        } else {
            istringstream iss(line);
            string token;
            vector<string> dataLine;

            while (getline(iss, token, ',')) {
                dataLine.push_back(token);
            }

            if (dataLine.empty()) {
                continue;
            }
            
            if (typeSpec == "nodes\r") {
                if (dataLine.size() != 3) {
                    cerr << "Warning. Invalid input for nodes." << endl;
                } else {
                    double x = stod(dataLine[0]);
                    double y = stod(dataLine[1]);
                    double z = stod(dataLine[2]);
                    Nodes->push_back(Vector3d(x, y, z));
                }
            } else if (typeSpec == "facenodes\r") {
                if (dataLine.size() != 3) {
                    cerr << "Warning. Invalid input for edges." << endl;
                } else {
                    int node1 = stoi(dataLine[0]);
                    int node2 = stoi(dataLine[1]);
                    int node3 = stoi(dataLine[2]);
                    facenodes->push_back({node1, node2, node3});
                }
            }
        }
    }

    inputFile.close();
}


void geometry_for_shell(vector<Vector3d>* Nodes, vector<vector<int> >* Face_Nodes, \
    vector<Vector3d>* Edges, vector<vector<int> >* Face_Edges, vector<vector<int> >* sign_face_edges, \
    vector<vector<int> >* EdgeIsBet, vector<vector<int> >* HingeIsBet, \
    vector<Vector3d>* edge_avg_normals) {

    // Input:
    // Nodes, Face_nodes
    //
    // Output:
    // Edges, Face_Edges, sign_face_edges,EdgeIsBet, HingeIsBet, edge_avg_normals

    int n_nodes = Nodes->size();
    int n_faces = Face_Nodes->size();

    int edge_index = 0;
    int hinge_index = 0;

    vector<int> third_node;

    for (int c = 0; c < n_faces; ++c) {
        int node1_number = (*Face_Nodes)[c][0];
        int node2_number = (*Face_Nodes)[c][1];
        int node3_number = (*Face_Nodes)[c][2];

        Vector3d node1_position = (*Nodes)[node1_number];
        Vector3d node2_position = (*Nodes)[node2_number];
        Vector3d node3_position = (*Nodes)[node3_number];

        Vector3d edge1 = node2_position - node1_position;
        Vector3d edge2 = node3_position - node2_position;
        Vector3d edge3 = node1_position - node3_position;

        Vector3d face_normal = edge3.cross(edge1);

        Vector3d face_unit_norm = face_normal.normalized();

        vector<int> edge1_between{ node1_number, node2_number };
        vector<int> edge1_between_negative{ node2_number, node1_number };

        vector<int> edge2_between{ node2_number, node3_number };
        vector<int> edge2_between_negative{ node3_number, node2_number };

        vector<int> edge3_between{ node3_number, node1_number };
        vector<int> edge3_between_negative{ node1_number, node3_number };

        //////////////////////////////////////////////////////////////
        /* Edge1 */
        bool bool_edge1_between_absent = true;
        bool bool_edge1_between_negative_absent = true;
        //int edge1_between_present = 0;
        int Locb1_between = -1;

        // check if edge1 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet->size(); i++) {
            if (edge1_between == (*EdgeIsBet)[i]) {
                bool_edge1_between_absent = false;
                Locb1_between = i;
            }
            else if (edge1_between_negative == (*EdgeIsBet)[i]) {
                bool_edge1_between_negative_absent = false;
                Locb1_between = i;
            }
        }
        // if not present, it is just an edge for now
        if (bool_edge1_between_absent && bool_edge1_between_negative_absent) {
            Edges->push_back(edge1);
            EdgeIsBet->push_back(edge1_between);
            third_node.push_back(node3_number);

            (*Face_Edges)[c].push_back(edge_index);
            edge_avg_normals->push_back(face_unit_norm);
            (*sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb1_between];
            int third_node_new = node3_number;

            HingeIsBet->push_back({ node1_number, node2_number, third_node_old, third_node_new });
            (*Face_Edges)[c].push_back(Locb1_between);
            edge_avg_normals->push_back(((*edge_avg_normals)[Locb1_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge1_between_absent) {
                (*sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge1_between_negative_absent) {
                (*sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }

        }
        //////////////////////////////////////////////////////////////
        /* Edge2 */
        bool bool_edge2_between_absent = true;
        bool bool_edge2_between_negative_absent = true;
        int Locb2_between = -1;

        // check if edge2 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet->size(); i++) {
            if (edge2_between == (*EdgeIsBet)[i]) {
                bool_edge2_between_absent = false;
                Locb2_between = i;
            }
            else if (edge2_between_negative == (*EdgeIsBet)[i]) {
                bool_edge2_between_negative_absent = false;
                Locb2_between = i;
            }
        }

        // if not present, it is just an edge for now
        if (bool_edge2_between_absent && bool_edge2_between_negative_absent) {
            Edges->push_back(edge2);
            EdgeIsBet->push_back(edge2_between);
            third_node.push_back(node1_number);

            (*Face_Edges)[c].push_back(edge_index);
            edge_avg_normals->push_back(face_unit_norm);
            (*sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb2_between];
            int third_node_new = node1_number;

            HingeIsBet->push_back({ node2_number, node3_number, third_node_old, third_node_new });
            (*Face_Edges)[c].push_back(Locb2_between);
            edge_avg_normals->push_back(((*edge_avg_normals)[Locb2_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge2_between_absent) {
                (*sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge2_between_negative_absent) {
                (*sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }
        }
        //////////////////////////////////////////////////////////////
        /* Edge3 */
        bool bool_edge3_between_absent = true;
        bool bool_edge3_between_negative_absent = true;
        int Locb3_between = -1;

        // check if edge2 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet->size(); i++) {
            if (edge3_between == (*EdgeIsBet)[i]) {
                bool_edge3_between_absent = false;
                Locb3_between = i;
            }
            else if (edge3_between_negative == (*EdgeIsBet)[i]) {
                bool_edge3_between_negative_absent = false;
                Locb3_between = i;
            }
        }

        // if not present, it is just an edge for now
        if (bool_edge3_between_absent && bool_edge3_between_negative_absent) {
            Edges->push_back(edge3);
            EdgeIsBet->push_back(edge3_between);
            third_node.push_back(node2_number);

            (*Face_Edges)[c].push_back(edge_index);
            edge_avg_normals->push_back(face_unit_norm);
            (*sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb3_between];
            int third_node_new = node2_number;

            HingeIsBet->push_back({ node3_number, node1_number, third_node_old, third_node_new });
            (*Face_Edges)[c].push_back(Locb3_between);
            edge_avg_normals->push_back(((*edge_avg_normals)[Locb3_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge3_between_absent) {
                (*sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge3_between_negative_absent) {
                (*sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }
        }
    }

    int n_edges = edge_index;
    int n_hinges = hinge_index;

    // check for debugging: check if size of edges and hinges vecs are correct
    if (EdgeIsBet->size() != n_edges || HingeIsBet->size() != n_hinges) {
        cout << "error in edge/hinges vectors" << endl;
    }
    return;
}

void Calculate_fi_ti_ci_shape_operator(Vector3d pi, Vector3d pj, Vector3d pk, double xi_i, double xi_j, double xi_k, \
    int s_i, int s_j, int s_k, Vector3d tau_i0, Vector3d tau_j0, Vector3d tau_k0, \
    vector<vector<double> >* init_fs, vector<vector<double> >* init_cs, vector<vector<Vector3d> >* init_ts, \
    vector<Matrix<double, 3, 3> >* lambdas) {

    Vector3d vi = pk - pj;
    Vector3d vj = pi - pk;
    Vector3d vk = pj - pi;

    // edge lengths
    double li = vi.norm();
    double lj = vj.norm();
    double lk = vk.norm();

    // triangle face normal
    Vector3d normal = vk.cross(vi);
    double A = normal.norm() / 2;
    Vector3d unit_norm = normal.normalized();

    // t_i's (tangent)
    Vector3d t_i = vi.cross(unit_norm);
    Vector3d t_j = vj.cross(unit_norm);;
    Vector3d t_k = vk.cross(unit_norm);;

    // c_i's : scalars
    double c_i = 1 / (A * li * t_i.dot(tau_i0) / t_i.norm());
    double c_j = 1 / (A * lj * t_j.dot(tau_j0) / t_j.norm());
    double c_k = 1 / (A * lk * t_k.dot(tau_k0) / t_k.norm());

    // f_i's : scalars
    double f_i = unit_norm.dot(tau_i0);
    double f_j = unit_norm.dot(tau_j0);
    double f_k = unit_norm.dot(tau_k0);

    // Combining into larger containers 

    int s[3] = { 0,0,0 };
    double xi[3] = { 0,0,0 };
    Vector3d tau_0[3];
    vector<double> f(3, 0);
    vector<double> c(3, 0);

    Vector3d temp;
    temp << 0, 0, 0;
    vector<Vector3d> t(3, temp);

    f[0] = f_i;
    f[1] = f_j;
    f[2] = f_k;

    c[0] = c_i;
    c[1] = c_j;
    c[2] = c_k;

    t[0] = t_i;
    t[1] = t_j;
    t[2] = t_k;

    xi[0] = xi_i;
    xi[1] = xi_j;
    xi[2] = xi_k;

    s[0] = s_i;
    s[1] = s_j;
    s[2] = s_k;

    tau_0[0] = tau_i0;
    tau_0[1] = tau_j0;
    tau_0[2] = tau_k0;

    MatrixXd lambda = MatrixXd::Zero(3, 3);

    // Shape operator calculation
    for (int i = 0; i < 3; i++) {
        lambda += (s[i] * xi[i] - f[i]) * c[i] * (t[i] * (t[i].transpose()));
    }

    // pushing into the larger vectors
    init_fs->push_back(f);
    init_cs->push_back(c);
    init_ts->push_back(t);
    lambdas->push_back(lambda);

}

void initial_geometric_params(vector<Vector3d>* Nodes, vector<Vector3d>* Edges, vector<vector<int> >* face_nodes, \
    vector<vector<int> >* face_edges, vector<vector<int> >* sign_face_edges, vector<Vector3d>* edge_n_avg, \
    vector<double>* q, vector<double>* undef_el, vector<Vector3d>* tau_0, vector<double>* xi_s, \
    vector<vector<double> >* init_fs, vector<vector<double> >* init_cs, vector<vector<Vector3d> >* init_ts, \
    vector<Matrix<double, 3, 3> >* shape_operators_faces) {
    // calculating undeformed edge lengths, tau_0s for each edge, initial xi's, intitial dof vector and 
    // for initial curvature: initial t's, intial f's, initial c's
    // shape operator 3*3 matrices for each face
    // intial DoF vector

    int n_edges = Edges->size();
    int n_nodes = Nodes->size();
    int n_faces = face_nodes->size();

    for (int i = 0; i < n_nodes; i++) {

        (*q).push_back((*Nodes)[i](0));
        (*q).push_back((*Nodes)[i](1));
        (*q).push_back((*Nodes)[i](2));
    }

    for (int i = 0; i < n_edges; i++) {
        q->push_back(0);
        undef_el->push_back((*Edges)[i].norm());
        tau_0->push_back((*Edges)[i].cross((*edge_n_avg)[i]));
        xi_s->push_back(0); // 0 intially
    }

    int Face_i_nodes[3] = { -1,-1,-1 };
    int Face_i_edges[3] = { -1,-1,-1 };
    int s_is[3] = { 0, 0, 0 };
    Vector3d p_is[3];
    Vector3d tau_0_is[3];
    double xi_is[3] = { 0,0,0 };

    for (int i = 0; i < n_faces; i++) {
        for (int j = 0; j < 3; j++) {
            Face_i_nodes[j] = (face_nodes->at(i))[j];
            Face_i_edges[j] = (*face_edges)[i][j];
        }

        for (int j = 0; j < 3; j++) {
            p_is[j] << (*q)[(3 * Face_i_nodes[j])], (*q)[(3 * Face_i_nodes[j]) + 1], (*q)[(3 * Face_i_nodes[j]) + 2];
            xi_is[j] = (*q)[(3 * n_nodes) + Face_i_edges[j]];
            tau_0_is[j] << (*tau_0)[Face_i_edges[j]];
            s_is[j] = (*sign_face_edges)[i][j];
        }

        Calculate_fi_ti_ci_shape_operator(p_is[0], p_is[1], p_is[2], xi_is[0], xi_is[1], xi_is[2], s_is[0], s_is[1], s_is[2], \
            tau_0_is[0], tau_0_is[1], tau_0_is[2], init_fs, init_cs, init_ts, shape_operators_faces);
    }

}

void MassMatrix_and_GravityForce(const double totalM, const double totalL, const double rho, int n_nodes, int n_edges, int n_dof, const double g[3], \
    MatrixXd* massMat, VectorXd* F_g) {
    // calculate massMatrix(n_dof*n_dof) and gravitational Force vector (n_dof*1)

    double unit_mass_node = 1; // can be calculated using material properties
    vector<double> massVec;

    for (int i = 0; i < 3 * n_nodes; i++) {
        massVec.push_back(unit_mass_node);
    }
    for (int j = 0; j < n_edges; j++) {
        massVec.push_back(0);
    }
    VectorXd massVecXd(n_dof);

    for (int k = 0; k < n_dof; k++) {
        massVecXd[k] = massVec[k];
    }
    (*massMat) = massVecXd.asDiagonal();

    VectorXd F_g_temp(n_dof);
    for (int i = 0; i < n_nodes; i++) {
        F_g_temp[3 * i] = massVec[3 * i] * g[0];
        F_g_temp[3 * i + 1] = massVec[3 * i + 1] * g[1];
        F_g_temp[3 * i + 2] = massVec[3 * i + 2] * g[2];
    }
    (*F_g) = F_g_temp;
}

void SetBoundaryConditions(vector<int> Fixed_node_indices, vector<int> Fixed_edge_indices, int n_nodes, int n_dof, \
    vector<int>* fixedDOF, vector<int>* freeDOF) {
    // set fixed and free DoF from user input of fixed nodes and edges
    for (int i = 0; i < Fixed_node_indices.size(); i++) {
        for (int j = 0; j < 3; j++) {
            (*fixedDOF).push_back(3 * Fixed_node_indices[i] + j);
        }
    }

    for (int k = 0; k < Fixed_edge_indices.size(); k++) {
        (*fixedDOF).push_back(n_nodes * 3 + Fixed_edge_indices[k]);
    }

    vector<int> dummy(n_dof, 1);
    for (int i : (*fixedDOF)) {
        dummy[i] = 0;
    }

    for (int i = 0; i < n_dof; ++i) {
        if (dummy[i] == 1) {
            (*freeDOF).push_back(i);
        }
    }
}

/* Elastic Bending Force calculation functions */

RowVector3d delfi_by_delpk(Vector3d tau_i0, Vector3d t_k, Vector3d unit_norm, double A) {
    // function needed for bending energy gradient
    Vector3d temp = tau_i0.dot(t_k) * unit_norm / (2 * A);
    RowVector3d del_fi_del_pk = temp.transpose();
    return del_fi_del_pk;
}

Matrix<double, 3, 3> ddel_fi_by_del_pk1_pk2(Vector3d vi, Vector3d vj, Vector3d vk, Vector3d tau_i0, \
    Vector3d unit_norm, double A, char k1, char k2) {
    // function needed for bending energy hessian
    Vector3d vk1, vk2;
    if (k1 == 'i') vk1 = vi;
    else if (k1 == 'j') vk1 = vj;
    else if (k1 == 'k') vk1 = vk;
    else cout << "error: k1 should be either i, j, or k";

    Vector3d tk1 = vk1.cross(unit_norm);

    if (k2 == 'i') vk2 = vi;
    else if (k2 == 'j') vk2 = vj;
    else if (k2 == 'k') vk2 = vk;
    else cout << "error: k2 should be either i, j, or k";

    Vector3d tk2 = vk2.cross(unit_norm);

    Matrix<double, 3, 3> ddel_fi_pk1_pk2;
    ddel_fi_pk1_pk2 = (1 / (4 * pow(A, 2))) * (tau_i0.dot(tk1) * ((unit_norm * (tk2.transpose())) + (tk2 * (unit_norm.transpose()))));
    return ddel_fi_pk1_pk2;
}



void grad_hess_Eb_shell_midedgeNormal(double Kbend, Vector3d pi, Vector3d pj, Vector3d pk, double xi_i, double xi_j, double xi_k, \
    int s_i, int s_j, int s_k, Vector3d tau_i0, Vector3d tau_j0, Vector3d tau_k0, \
    vector<Vector3d> init_ts, vector<double> init_cs, vector<double> init_fs, \
    Vector<double, 12>& dF, Matrix<double, 12, 12>& dJ)
{
    Vector3d vi = pk - pj;
    Vector3d vj = pi - pk;
    Vector3d vk = pj - pi;

    // edge lengths
    double li = vi.norm();
    double lj = vj.norm();
    double lk = vk.norm();

    // triangle face normal
    Vector3d normal = vk.cross(vi);
    double A = normal.norm() / 2;
    Vector3d unit_norm = normal.normalized();

    // t_i's (tangent)
    Vector3d t_i = vi.cross(unit_norm);
    Vector3d t_j = vj.cross(unit_norm);;
    Vector3d t_k = vk.cross(unit_norm);;

    // c_i's : scalars
    double c_i = 1 / (A * li * t_i.dot(tau_i0) / t_i.norm());
    double c_j = 1 / (A * lj * t_j.dot(tau_j0) / t_j.norm());
    double c_k = 1 / (A * lk * t_k.dot(tau_k0) / t_k.norm());

    // f_i's : scalars
    double f_i = unit_norm.dot(tau_i0);
    double f_j = unit_norm.dot(tau_j0);
    double f_k = unit_norm.dot(tau_k0);

    // Combining into larger containers 
    double f[3], c[3], xi[3];
    int s[3];
    Vector3d t[3], tau_0[3];

    f[0] = f_i;
    f[1] = f_j;
    f[2] = f_k;

    c[0] = c_i;
    c[1] = c_j;
    c[2] = c_k;

    xi[0] = xi_i;
    xi[1] = xi_j;
    xi[2] = xi_k;

    s[0] = s_i;
    s[1] = s_j;
    s[2] = s_k;

    t[0] = t_i;
    t[1] = t_j;
    t[2] = t_k;

    tau_0[0] = tau_i0;
    tau_0[1] = tau_j0;
    tau_0[2] = tau_k0;

    ///////////////// Gradient of Energy

    // initialize
    RowVector3d del_E_del_pi = RowVector3d::Zero();
    RowVector3d del_E_del_pj = RowVector3d::Zero();
    RowVector3d del_E_del_pk = RowVector3d::Zero();

    double del_E_del_xi_i = 0;
    double del_E_del_xi_j = 0;
    double del_E_del_xi_k = 0;

    // derivatives wrt pi's
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            del_E_del_pi = del_E_del_pi - c[i] * c[j] * ((s[i] * xi[i] - f[i]) * delfi_by_delpk(tau_0[j], t_i, unit_norm, A) + \
                (s[j] * xi[j] - f[j]) * delfi_by_delpk(tau_0[i], t_i, unit_norm, A)) * pow((t[i].dot(t[j])), 2) - \
                2 * c[i] * init_cs[j] * init_fs[j] * delfi_by_delpk(tau_0[i], t_i, unit_norm, A) * pow((t[i].dot(init_ts[j])), 2);


            del_E_del_pj = del_E_del_pj - c[i] * c[j] * ((s[i] * xi[i] - f[i]) * delfi_by_delpk(tau_0[j], t_j, unit_norm, A) + \
                (s[j] * xi[j] - f[j]) * delfi_by_delpk(tau_0[i], t_j, unit_norm, A)) * pow((t[i].dot(t[j])), 2) - \
                2 * c[i] * init_cs[j] * init_fs[j] * delfi_by_delpk(tau_0[i], t_j, unit_norm, A) * pow((t[i].dot(init_ts[j])), 2);


            del_E_del_pk = del_E_del_pk - c[i] * c[j] * ((s[i] * xi[i] - f[i]) * delfi_by_delpk(tau_0[j], t_k, unit_norm, A) + \
                (s[j] * xi[j] - f[j]) * delfi_by_delpk(tau_0[i], t_k, unit_norm, A)) * pow((t[i].dot(t[j])), 2) - \
                2 * c[i] * init_cs[j] * init_fs[j] * delfi_by_delpk(tau_0[i], t_k, unit_norm, A) * pow((t[i].dot(init_ts[j])), 2);
        }
    }
    // derivatives wrt xi_i's

    for (int j = 0; j < 3; j++) {
        del_E_del_xi_i += (2 * c_i * s_i) * c[j] * (s[j] * xi[j] - f[j]) * pow(t_i.dot(t[j]), 2) + \
            2 * c_i * s_i * init_cs[j] * init_fs[j] * pow(t_i.dot(init_ts[j]), 2);

        del_E_del_xi_j += (2 * c_j * s_j) * c[j] * (s[j] * xi[j] - f[j]) * pow(t_j.dot(t[j]), 2) + \
            2 * c_j * s_j * init_cs[j] * init_fs[j] * pow(t_j.dot(init_ts[j]), 2);

        del_E_del_xi_k += (2 * c_k * s_k) * c[j] * (s[j] * xi[j] - f[j]) * pow(t_k.dot(t[j]), 2) + \
            2 * c_k * s_k * init_cs[j] * init_fs[j] * pow(t_k.dot(init_ts[j]), 2);
    }

    // collect in the gradient vector in the sequence: del_by [pi, pj, pk, xi_i, xi_j, xi_k]   size = (1*12)

    RowVectorXd grad(12);

    grad << del_E_del_pi, del_E_del_pj, del_E_del_pk, del_E_del_xi_i, del_E_del_xi_j, del_E_del_xi_k;
    // dF = Kbend * grad;
     dF = Kbend * grad.transpose();
    //////////////// Hessian of Energy
    // 
    // ddel_E_by_del_xi1_xi2
    MatrixXd ddel_E_del_xis = MatrixXd::Zero(3, 3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ddel_E_del_xis(i, j) = 2 * (c[i] * c[j]) * (s[i] * s[j]) * pow(t[i].dot(t[j]), 2);
        }
    }
    // ddel_E_by_del_p1_xi2
    MatrixXd ddel_E_del_xi_del_p_i = MatrixXd::Zero(3, 3);
    MatrixXd ddel_E_del_xi_del_p_j = MatrixXd::Zero(3, 3);
    MatrixXd ddel_E_del_xi_del_p_k = MatrixXd::Zero(3, 3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            ddel_E_del_xi_del_p_i.row(i) += -(2 * (c[i] * s[i])) * c[j] * \
                delfi_by_delpk(tau_0[j], t_i, unit_norm, A) * pow(t[i].dot(t[j]), 2);

            ddel_E_del_xi_del_p_j.row(i) += -(2 * (c[i] * s[i])) * c[j] * \
                delfi_by_delpk(tau_0[j], t_j, unit_norm, A) * pow(t[i].dot(t[j]), 2);

            ddel_E_del_xi_del_p_k.row(i) += -(2 * (c[i] * s[i])) * c[j] * \
                delfi_by_delpk(tau_0[j], t_k, unit_norm, A) * pow(t[i].dot(t[j]), 2);

        }
    }

    Matrix<double, 3, 9> hess_xi_pi;
    hess_xi_pi << ddel_E_del_xi_del_p_i, ddel_E_del_xi_del_p_j, ddel_E_del_xi_del_p_k;

    // ddel_E_by_del_p1_p2
    //MatrixXd Id3 = MatrixXd::Identity(3, 3);
    Matrix<double, 3, 3> init_mat;
    init_mat = Matrix<double, 3, 3>::Zero();

    vector < Matrix<double, 3, 3> > ddel_E_del_ps(9, init_mat);

    char char_k1, char_k2;
    int num_k1, num_k2;

    for (int k = 0; k < 9; k++) {

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {

                if (k == 0) {
                    char_k1 = 'i';
                    char_k2 = 'i';
                    num_k1 = 0;
                    num_k2 = 0;
                }

                else if (k == 1) {
                    char_k1 = 'i';
                    char_k2 = 'j';
                    num_k1 = 0;
                    num_k2 = 1;
                }

                else if (k == 2) {
                    char_k1 = 'i';
                    char_k2 = 'k';
                    num_k1 = 0;
                    num_k2 = 2;
                }

                else if (k == 3) {
                    char_k1 = 'j';
                    char_k2 = 'i';
                    num_k1 = 1;
                    num_k2 = 0;
                }

                else if (k == 4) {
                    char_k1 = 'j';
                    char_k2 = 'j';
                    num_k1 = 1;
                    num_k2 = 1;
                }

                else if (k == 5) {
                    char_k1 = 'j';
                    char_k2 = 'k';
                    num_k1 = 1;
                    num_k2 = 2;
                }

                else if (k == 6) {
                    char_k1 = 'k';
                    char_k2 = 'i';
                    num_k1 = 2;
                    num_k2 = 0;
                }

                else if (k == 7) {
                    char_k1 = 'k';
                    char_k2 = 'j';
                    num_k1 = 2;
                    num_k2 = 1;
                }

                else if (k == 8) {
                    char_k1 = 'k';
                    char_k2 = 'k';
                    num_k1 = 2;
                    num_k2 = 2;
                }

                else {
                    cout << "error in k: should be in {0,8}";
                }

                ddel_E_del_ps[k] += -(c[i] * c[j]) * pow(t[i].dot(t[j]), 2) * (\
                    (s[i] * xi[i] - f[i]) * ddel_fi_by_del_pk1_pk2(vi, vj, vk, tau_0[j], unit_norm, A, char_k2, char_k1) - \
                    (delfi_by_delpk(tau_0[j], t[num_k1], unit_norm, A).transpose() * delfi_by_delpk(tau_0[i], t[num_k2], unit_norm, A)) + \
                    (s[j] * xi[j] - f[j]) * ddel_fi_by_del_pk1_pk2(vi, vj, vk, tau_0[i], unit_norm, A, char_k2, char_k1) - \
                    (delfi_by_delpk(tau_0[i], t[num_k1], unit_norm, A).transpose() * delfi_by_delpk(tau_0[j], t[num_k2], unit_norm, A)))
                    + \
                    2 * c[i] * init_cs[j] * -init_fs[j] * ddel_fi_by_del_pk1_pk2(vi, vj, vk, tau_0[i], unit_norm, A, char_k2, char_k1) * pow(t[i].dot(init_ts[j]), 2);

            }
        }
    }
    // combining all into the hessian matrix (12*12)
    MatrixXd hessE(12, 12);
    hessE.block <9, 9>(0, 0) << ddel_E_del_ps[0], ddel_E_del_ps[1], ddel_E_del_ps[2], ddel_E_del_ps[3], ddel_E_del_ps[4], \
        ddel_E_del_ps[5], ddel_E_del_ps[6], ddel_E_del_ps[7], ddel_E_del_ps[8];
    hessE.block <9, 3>(0, 9) << hess_xi_pi.transpose();
    hessE.block <3, 9>(9, 0) << hess_xi_pi;
    hessE.block <3, 3>(9, 9) << ddel_E_del_xis;

     dJ = Kbend * hessE.transpose(); //correct
    // dJ = Kbend * hessE; // not correct
}

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
void get_bendingForce(vector<double>* q_pointer, int n_nodes, int n_faces, \
    vector <vector<int> >* face_nodes, vector<vector<int> >* face_edges, double Kbend, \
    vector <Vector3d>* tau_0, vector<vector<int> >* sign_faces, \
    vector< vector <Vector3d> >* init_ts, vector< vector<double> >* init_cs, vector< vector<double> >* init_fs, \
    VectorXd* Fb_pointer, MatrixXd* Jb_pointer) {

    int Face_i_nodes[3] = { -1,-1,-1 };
    int Face_i_edges[3] = { -1,-1,-1 };
    int s_is[3] = { 0, 0, 0 };
    Vector3d p_is[3];
    Vector3d tau_0_is[3];
    double xi_is[3] = { 0,0,0 };

    int n_dof = q_pointer->size();
    VectorXd Fb = VectorXd::Zero(n_dof);
    MatrixXd Jb = MatrixXd::Zero(n_dof, n_dof);

    for (int i = 0; i < n_faces; i++) {
        for (int j = 0; j < 3; j++) {
            Face_i_nodes[j] = (*face_nodes)[i][j];
            // cout<<Face_i_nodes[j]<<" ";
            Face_i_edges[j] = (*face_edges)[i][j];
            // cout<<Face_i_edges[j]<<" ";
        }

        for (int j = 0; j < 3; j++) {
            p_is[j] << (*q_pointer)[(3 * Face_i_nodes[j])], (*q_pointer)[(3 * Face_i_nodes[j]) + 1], (*q_pointer)[(3 * Face_i_nodes[j]) + 2];
            xi_is[j] = (*q_pointer)[(3 * n_nodes) + Face_i_edges[j]];
            tau_0_is[j] << (*tau_0)[Face_i_edges[j]];
            s_is[j] = (*sign_faces)[i][j];
        }

        Vector<double, 12> dFb = Vector<double,12>::Zero();
        Matrix<double, 12, 12> dJb = Matrix<double,12,12>::Zero();

        grad_hess_Eb_shell_midedgeNormal(Kbend, p_is[0], p_is[1], p_is[2], xi_is[0], xi_is[1], xi_is[2], s_is[0], s_is[1], s_is[2], \
            tau_0_is[0], tau_0_is[1], tau_0_is[2], (*init_ts)[i], (*init_cs)[i], (*init_fs)[i], dFb, dJb);

        int ind[12];
        for (int j = 0; j < 3; j++) {
            ind[3 * j] = (3 * Face_i_nodes[j]);
            ind[3 * j + 1] = (3 * Face_i_nodes[j] + 1);
            ind[3 * j + 2] = (3 * Face_i_nodes[j] + 2);
        }
        for (int j = 0; j < 3; j++) {
            ind[9 + j] = 3 * n_nodes + Face_i_edges[j];
        }
        // how to add in a force vector and jacobian matrix ndof*ndof
                // Fs(ind) -= dFs; // incorrect as not in sequence?
                // Js.block<12, 12>(ind[0], ind[0]) -= dJs; // incorrect as not in sequence
        int j = 0;
        for (int i : ind) {
            Fb(i) = Fb(i) - dFb(j);
            j++;
        }

        int p_counter = 0, q_counter = 0;
        for (int i : ind) {
            q_counter = 0;
            for (int j : ind) {
                Jb(i, j) =  Jb(i,j) - dJb(p_counter, q_counter);
                q_counter++;
            }
            p_counter++;
        }
        *Fb_pointer = Fb;
        *Jb_pointer = Jb;
    }
}


/* Elastic Stretching Force calculation functions */

void grad_hess_Es_shell(Vector3d node0, Vector3d node1, double undef_el, double Kstretch, \
    Vector<double,6>& dF, Matrix<double,6,6>& dJ)
{
    Vector3d edge = node1 - node0;
    double edgeLen = edge.norm();

    // gradient
    Vector3d tangent = edge / edgeLen;
    double epsX = edgeLen / undef_el - 1;
    Vector3d dF_unit = Kstretch / undef_el * tangent * epsX;

    //VectorXd dF(6);
    dF << -dF_unit, dF_unit;

    // hessian
    MatrixXd Id3 = Matrix<double, 3, 3>::Identity();
    MatrixXd M_temp(3, 3);
    M_temp = (Kstretch / undef_el) * ((1 / undef_el - 1 / edgeLen) * Id3 + \
        1 / edgeLen * (edge * edge.transpose()) / pow(edgeLen, 2));

    //MatrixXd dJ(6, 6);
    dJ << M_temp, -M_temp, -M_temp, M_temp;
}


void get_stretchingForce(vector<double>* q_pointer, int n_edges, vector<vector<int> >* EdgeIsBetn_pointer, \
    vector<double>* undef_els, double Kstretch, VectorXd* Fs_pointer, MatrixXd* Js_pointer)
{
    int n_dof = q_pointer->size();
    VectorXd Fs = VectorXd::Zero(n_dof);
    MatrixXd Js = MatrixXd::Zero(n_dof, n_dof);

    for (int i = 0; i < n_edges; i++) {
        int node0ind = (*EdgeIsBetn_pointer)[i][0];
        int node1ind = (*EdgeIsBetn_pointer)[i][1];

        Vector3d node0, node1;
        node0 << (*q_pointer)[3 * node0ind], (*q_pointer)[3 * node0ind + 1], (*q_pointer)[3 * node0ind + 2];
        node1 << (*q_pointer)[3 * node1ind], (*q_pointer)[3 * node1ind + 1], (*q_pointer)[3 * node1ind + 2];

        int ind[6] = { 3 * node0ind, 3 * node0ind + 1, 3 * node0ind + 2, \
            3 * node1ind, 3 * node1ind + 1, 3 * node1ind + 2 };

        Vector<double, 6> dFs = Vector<double,6>::Zero();
        Matrix<double, 6, 6> dJs = Matrix<double,6,6>::Zero();

        grad_hess_Es_shell(node0, node1, (*undef_els)[i], Kstretch, dFs, dJs);

        Fs(ind) -= dFs;
        int p_counter = 0, q_counter = 0;
        for (int i : ind) {
            q_counter = 0;
            for (int j : ind) {
                Js(i, j) =  Js(i,j) - dJs(p_counter, q_counter);
                q_counter++;
            }
            p_counter++;
        }
        // Js.block<6, 6>(ind[0], ind[0]) -= dJs; // this is incorrect
    }
    *Fs_pointer = Fs;
    *Js_pointer = Js;
}


// .........................................................................................
// main function

int main()
{
    vector<Vector3d>* Nodes = new vector<Vector3d>(0);
    vector<vector<int> >* Face_Nodes = new vector<vector<int> >(0);
    // Read the input txt file into:
    // Nodes
    // Face_Nodes

    string inputFileName = "two_triangle_input.txt"; // two triangles
    // string inputFileName = "input_shell_plate_cantilever_for_c++.txt"; // two triangles
    inputProcessorFunc(inputFileName, Nodes, Face_Nodes);

    // cout << "Node positions:" << endl;
    // for (int i = 0; i < Nodes->size(); i++) {
    //     cout << (*Nodes)[i] << endl;
    // }
    // cout << "Face Nodes:" << endl;
    // for (int i = 0; i < Face_Nodes->size(); i++) {
    //     for (int j = 0; j < (*Face_Nodes)[i].size(); j++) {
    //         cout << (*Face_Nodes)[i][j] << ", ";
    //     }
    //     cout << ";" << endl;
    // }
// return 0;
    // some required numbers:
    int n_nodes = Nodes->size();
    int n_faces = Face_Nodes->size();

    // create initial structs:

    vector<Vector3d>* Edges = new vector<Vector3d>(0);
    vector<vector<int> >* Face_Edges = new vector<vector<int> >(n_faces);
    vector<vector<int> >* sign_face_edges = new vector<vector<int> >(n_faces);
    vector<vector<int> >* EdgeIsBet = new vector<vector<int> >(0);
    vector<vector<int> >* HingeIsBet = new vector<vector<int> >(0);
    vector<Vector3d>* edge_avg_normals = new vector<Vector3d>(0);

    geometry_for_shell(Nodes, Face_Nodes, \
        Edges, Face_Edges, sign_face_edges, EdgeIsBet, HingeIsBet, edge_avg_normals);

    // cout << "Edge vectors:" << endl;
    // for (int i = 0; i < Edges->size(); i++) {
    //     cout << (*Edges)[i] << endl;
    // }
    // cout << "Face Edges:" << endl;
    // for (int i = 0; i < Face_Edges->size(); i++) {
    //     for (int j = 0; j < (*Face_Edges)[i].size(); j++) {
    //         cout << (*Face_Edges)[i][j] << " " ;
    //     }
    //     cout << endl;
    // }
    // cout << "Sign face edges:" << endl;
    // for (int i = 0; i < sign_face_edges->size(); i++) {
    //     for (int j = 0; j < (*sign_face_edges)[i].size(); j++) {
    //         cout << (*sign_face_edges)[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "Edge avg normals:" << endl;
    // for (int i = 0; i < edge_avg_normals->size(); i++) {
    //     cout << (*edge_avg_normals)[i] << endl;
    // }


    // create other structs:
    vector<double>* q = new vector<double>(0);
    vector<double>* undef_el = new vector<double>(0);
    vector<Vector3d>* tau_0 = new vector<Vector3d>(0);
    vector<double>* xi_s = new vector<double>(0);
    vector<vector<double> >* init_fs = new vector<vector<double> >(0);
    vector<vector<double> >* init_cs = new vector<vector<double> >(0);
    vector<vector<Vector3d> >* init_ts = new vector<vector<Vector3d> >(0);
    vector<Matrix<double, 3, 3> >* shape_operators_faces = new vector<Matrix<double, 3, 3> >(0);

    initial_geometric_params(Nodes, Edges, Face_Nodes, \
        Face_Edges, sign_face_edges, edge_avg_normals, \
        q, undef_el, tau_0, xi_s, \
        init_fs, init_cs, init_ts, \
        shape_operators_faces);

     // some required numbers
    int n_dof = q->size();
    int n_edges = Edges->size();

    // cout << "q intially:" << endl;
    // for (int i = 0; i < q->size(); i++) {
    //     cout << (*q)[i] << endl;
    // }
    

    // material properties
    const double Y = pow(10, 6); // N/m^2 (shell Young's Modulus)
    const double rho = 1200; // kg/m^3 (shell material density)
    const double nu = 0.5; // (poisson ratio)
    const double thickness = 0.001; // m (shell thickness)
    const double edge_len = 0.01; //m (avg edge length)
    const double totalM = 1; // kg (need to correct)
    const double totalL = 1; //m (need to correct)
    const double g[3] = { 0,0,9.8 };

    double Kstretch = 0.5 * sqrt(3) * Y * thickness * pow(edge_len, 2);
    double Kbend = 2 * Y * pow(thickness, 3) / (sqrt(3) * 12);
    // for now
    Kstretch = 1000;
    Kbend = 1;

    // creating other structs

    VectorXd Fg = VectorXd::Zero(n_dof) ;
    MatrixXd massMat = MatrixXd::Zero(n_dof,n_dof);

    MassMatrix_and_GravityForce(totalM, totalL, rho, n_nodes, n_edges, n_dof, g, \
        & massMat, &Fg);

    // Fixed and free dofs
    vector<int>* fixedDOF = new vector<int>(0);
    vector<int>* freeDOF = new vector<int>(0);
    vector<int> Fixed_node_indices{ 0,1,2 };
    vector<int> Fixed_edge_indices{ 0,1,2 };

    SetBoundaryConditions(Fixed_node_indices, Fixed_edge_indices, n_nodes, n_dof, \
        fixedDOF, freeDOF);

    int n_freeDOF = freeDOF->size();

    /* Simulation */

    // Simulation parameters
    double dt = 0.001; // sec
    double totalTime = 0.1; // sec
    double ctime = 0; //current time (sec)

    int Nsteps = round(totalTime / dt);

    VectorXd q_vectorXd = VectorXd::Zero(n_dof);
    VectorXd q0_vectorXd = VectorXd::Zero(n_dof);
    VectorXd u_vectorXd = VectorXd::Zero(n_dof);

    vector<double>* q0 = new vector<double>(n_dof); // intial DOF vector
    *q0 = *q;

    vector<double>* u = new vector<double>(n_dof, 0.0);

    // creating storage structs 

    VectorXd F_bending = VectorXd::Zero(n_dof);
    VectorXd F_stretching = VectorXd::Zero(n_dof);
    VectorXd F_total = VectorXd::Zero(n_dof);

    MatrixXd J_bending = MatrixXd::Zero(n_dof, n_dof);
    MatrixXd J_stretching = MatrixXd::Zero(n_dof, n_dof);
    MatrixXd J_total = MatrixXd::Zero(n_dof, n_dof);

    VectorXd f = VectorXd::Zero(n_dof);
    MatrixXd J = MatrixXd::Zero(n_dof, n_dof);

    VectorXd f_free = VectorXd::Zero(n_freeDOF);
    MatrixXd J_free = MatrixXd::Zero(n_freeDOF, n_freeDOF);

    VectorXd dq_free = VectorXd::Zero(n_freeDOF);
    VectorXd q_free = VectorXd::Zero(n_freeDOF);

    int iter = 0;
    double tol = 0.01;
    double error = 10 * tol;

    // get_stretchingForce(q, n_edges, EdgeIsBet, undef_el, Kstretch, &F_stretching, &J_stretching);


    // get_bendingForce(q, n_nodes, n_faces, Face_Nodes, Face_Edges, Kbend, tau_0, sign_face_edges, \
    //     init_ts, init_cs, init_fs, &F_bending, &J_bending);


    // Forces and Jacobians initially
    // cout << "F_stretching is:" << endl;
    // cout << F_stretching << endl;
    // cout << "F_bending is:" << endl;
    // cout << F_bending << endl;

    // cout<< "J_stretching is:"<<endl;
    // cout<< J_stretching<<endl;
    // cout<< "J_bending is:"<<endl;
    // cout<< J_bending<<endl;



    for (int timeStep = 0; timeStep < Nsteps; timeStep++) {
        // current dof vector is guessed to be same as old
        *q = *q0;

        error = 10 * tol;

        iter = 0;

        while (error > tol) {

            get_stretchingForce(q, n_edges, EdgeIsBet, undef_el, Kstretch, &F_stretching, &J_stretching);


            get_bendingForce(q, n_nodes, n_faces, Face_Nodes, Face_Edges, Kbend, tau_0, sign_face_edges, \
                init_ts, init_cs, init_fs, &F_bending, &J_bending);


            F_total = F_bending + F_stretching + Fg;
            J_total = J_bending + J_stretching;

            // convert q to vector form for matrix vector multiplications
            for (int i = 0; i < n_dof; i++) {
                q_vectorXd(i) = (*q)[i];
                q0_vectorXd(i) = (*q0)[i];
                u_vectorXd(i) = (*u)[i];
            }

            // EOM
            f = ((1 / dt) * massMat) * ((1 / dt) * (q_vectorXd - q0_vectorXd) - u_vectorXd) - F_total;
            // Jacobian
            J = ((1 / pow(dt, 2)) * massMat) - J_total;

            for (int i = 0; i < freeDOF->size(); i++) {
                f_free(i) = f((*freeDOF)[i]);
                for (int j = 0; j < freeDOF->size(); j++) {
                    J_free(i, j) = J((*freeDOF)[i], (*freeDOF)[j]);
                }
            }
            // Newton Raphson Root finding method

            dq_free = J_free.inverse() * f_free; // dq_free is a n_freeDOF*1 VectorXd

            // update the vector<double>* q using dq_free
            for (int i = 0; i < freeDOF->size(); i++) {
                q->at((*freeDOF)[i]) -= dq_free(i);
            }
            
            error = f_free.lpNorm<1>();

            cout << "iter = " << iter << ", " << "error = " << error << endl;

            iter++;

        }

        // update velocity
        for (int i = 0; i < n_dof; i++) {
            (*u)[i] = (1 / dt) * ((*q)[i] - (*q0)[i]);
        }
        // store current dof vector into old
        *q0 = *q;

        ctime = ctime + dt;

        cout << "current time: " << ctime << endl;

    }

    /////////////////////////////////////////////////////////////////////

    return 0;

}
