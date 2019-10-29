#include "viewer.h"
#include <iomanip>

using namespace surface_mesh;
typedef Surface_mesh Mesh;

bool isClose(float a, float b, float tol) {
    return std::abs(a - b) <= tol;
}

bool isClose(const Point &a, const Point &b, float tol=1e-2f) {
    return isClose(a[0], b[0], tol) && isClose(a[1], b[1], tol) && isClose(a[2], b[2], tol);
}

void print(const Point &p, bool brackets=true) {
    if (brackets) {
        std::cout << "[" << p[0] << ", " << p[1] << ", " << p[2] << "]" << std::endl;
    } else {
        std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
    }
}

mesh_processing::MeshProcessing *createNewMeshProcessing() {
    mesh_processing::MeshProcessing* mp = new mesh_processing::MeshProcessing("../data/bunny.off", true);
    Surface_mesh &mesh = mp->mesh_;
    if (!(mesh.n_vertices() == 34835 && mesh.n_faces() == 69666)) {
        std::cout << "Wrong mesh was loaded." << std::endl;
        return nullptr;
    }
    return mp;
}

Mesh::Vertex interestingVertex(const Surface_mesh *mesh) {
    // Single out an interesting vertex
    Mesh::Vertex v_test;
    size_t k = 0;
    for (auto v: mesh->vertices()) {
        if (k == 101) {
            v_test = v;
        }
        k++;
    }
    return v_test;
}

int main(int argc, char **argv) {
    cout.precision(15);
    mesh_processing::MeshProcessing *mp = nullptr;
    Surface_mesh *mesh = nullptr;
    Mesh::Vertex v_test;
    Point p_test, p_ref;

    int n_successes = 0;
    int n_tests = 0;

    int n_successes_optional = 0;
    int n_tests_optional = 0;

    std::cout << "======================================================================" << std::endl;
    std::cout << "Smoothing test" << std::endl;
    std::cout << "======================================================================" << std::endl;


    std::cout << "Part 1) Explicit Smoothing" << std::endl;

    std::cout << "  1.1: Uniform weights ....................................... ";

    p_ref = {-0.0406255386769772, 0.0344971679151058, 0.0191147085279226};

    mp = createNewMeshProcessing();
    mesh = &mp->mesh_;
    v_test = interestingVertex(mesh);

    mp->uniform_smooth(10);
    p_test = mesh->position(v_test);

    if (isClose(p_test, p_ref, 1e-4f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(p_test);
    std::cout << "    reference: "; print(p_ref);
    n_tests++;

    delete mp;

    std::cout << "  1.2: Cotan weights ......................................... ";

    p_ref = {-0.036231093108654, 0.0348517224192619, 0.0214509889483452};

    mp = createNewMeshProcessing();
    mesh = &mp->mesh_;
    v_test = interestingVertex(mesh);

    mp->smooth(10);
    p_test = mesh->position(v_test);

    if (isClose(p_test, p_ref, 1e-4f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(p_test);
    std::cout << "    reference: "; print(p_ref);
    n_tests++;

    delete mp;

    std::cout << std::endl;
    std::cout << "Part 2) Implicit Smoothing ................................... ";

    p_ref = {-0.0364644825458527, 0.0358456671237946, 0.0209404285997152};

    mp = createNewMeshProcessing();
    mesh = &mp->mesh_;
    v_test = interestingVertex(mesh);

    mp->implicit_smoothing(1e-4f);
    p_test = mesh->position(v_test);

    if (isClose(p_test, p_ref, 1e-4f)) {
        std::cout << "PASS" << std::endl;
        n_successes++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(p_test);
    std::cout << "    reference: "; print(p_ref);
    n_tests++;

    delete mp;

    std::cout << std::endl;
    std::cout << "Part 3) Feature Enhancement (Optional)" << std::endl;

    std::cout << "  3.1: Uniform weights ....................................... ";

    p_ref = {-0.0315384604036808, 0.0357928313314915, 0.0238772910088301};

    mp = createNewMeshProcessing();
    mesh = &mp->mesh_;
    v_test = interestingVertex(mesh);

    mp->uniform_laplacian_enhance_feature(10, 2);
    p_test = mesh->position(v_test);

    if (isClose(p_test, p_ref, 1e-4f)) {
        std::cout << "PASS" << std::endl;
        n_successes_optional++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(p_test);
    std::cout << "    reference: "; print(p_ref);
    n_tests_optional++;

    delete mp;

    std::cout << "  3.2: Cotan weights ......................................... ";

    p_ref = {-0.0359329059720039, 0.0354382768273354, 0.0215410105884075};

    mp = createNewMeshProcessing();
    mesh = &mp->mesh_;
    v_test = interestingVertex(mesh);

    mp->cotan_laplacian_enhance_feature(10, 2);
    p_test = mesh->position(v_test);

    if (isClose(p_test, p_ref, 1e-4f)) {
        std::cout << "PASS" << std::endl;
        n_successes_optional++;
    } else {
        std::cout << "FAIL" << std::endl;
    }
    std::cout << "    value:     "; print(p_test);
    std::cout << "    reference: "; print(p_ref);
    n_tests_optional++;

    delete mp;

    std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "PASSED " << n_successes << "/" << n_tests << " TESTS" << std::endl;
    std::cout << "PASSED " << n_successes_optional << "/" << n_tests_optional << " OPTIONAL TESTS" << std::endl;
}
