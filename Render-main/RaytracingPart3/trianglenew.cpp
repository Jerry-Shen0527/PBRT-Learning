#include "trianglenew.h"
aabb Trianglenew::bounding_box() const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    /*const Point3f& p0 = mesh->p[v[0]];
    const Point3f& p1 = mesh->p[v[1]];
    const Point3f& p2 = mesh->p[v[2]];*/
    //cout << p[0].x << " " << p[0].y << endl;
    //cout << Union(aabb(p[0], p[1]), p[2]).maximum.z << endl;
    return Union(aabb(p[0], p[1]), p[2]);
}

bool Trianglenew::Intersect(const ray& ray, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    const Point3f& p0 = p[0];
    const Point3f& p1 = p[1];
    const Point3f& p2 = p[2];
    //cout << p[0].x << " " << p[0].y << endl;

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.origin());
    Point3f p1t = p1 - Vector3f(ray.origin());
    Point3f p2t = p2 - Vector3f(ray.origin());

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.direction()));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.direction(), kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);


    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < t_max * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det))
        return false;

    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
        (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
        std::abs(invDet);
    if (t <= deltaT) return false;

    // Compute triangle partial derivatives
    Vector3f dpdu, dpdv;
    Point2f uv[3];
    GetUVs(uv);

    // Compute deltas for triangle partial derivatives
    Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
    Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    bool degenerateUV = std::abs(determinant) < 1e-8;
    if (!degenerateUV) {
        Float invdet = 1 / determinant;
        dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
    if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
        // Handle zero determinant for triangle partial derivative matrix
        Vector3f ng = Cross(p2 - p0, p1 - p0);
        if (ng.LengthSquared() == 0)
            // The triangle is actually degenerate; the intersection is
            // bogus.
            return false;

        CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
    }

    // Compute error bounds for triangle intersection
    Float xAbsSum =
        (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
    Float yAbsSum =
        (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
    Float zAbsSum =
        (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
    Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

    // Interpolate $(u,v)$ parametric coordinates and hit point
    Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
    Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

    inter_record.p = pHit;

    Vector3f outward_normal(Normal3f(Normalize(Cross(dp02, dp12))));
    //Vector3f trans_N = ConvertVTrans(ObjectToWorld, outward_normal);
    inter_record.set_face_normal(ray, outward_normal);
    inter_record.mat_ptr = mat_ptr;
    inter_record.dpdu = dpdu;
    inter_record.dpdv = dpdv;
    inter_record.t = t;
    
}

bool Trianglenew::IntersectP(const ray& ray, double t_min, double t_max) const {
    const Point3f& p0 = p[0];
    const Point3f& p1 = p[1];
    const Point3f& p2 = p[2];

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.origin());
    Point3f p1t = p1 - Vector3f(ray.origin());
    Point3f p2t = p2 - Vector3f(ray.origin());

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.direction()));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.direction(), kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);


    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < t_max * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det))
        return false;

    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
        (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
        std::abs(invDet);
    if (t <= deltaT) return false;
    return true;
}

TriangleMeshNew::TriangleMeshNew(
    const Transform& ObjectToWorld, bool reverseOrientation, const int nTriangles,
    const int nVertices, const vector<Point3f>& v, vector<int>& f, vector< shared_ptr<material>>& m)
    : nTriangles(nTriangles),
    nVertices(nVertices),
    ObjectToWorld(ObjectToWorld),
    reverseOrientation(reverseOrientation)
{
    faceVID.assign(f.begin(), f.end());
    vPos.assign(v.begin(), v.end());
    mat_ptr.assign(m.begin(), m.end());

}

std::vector<std::shared_ptr<Primitive>> CreateTriangleMesh(
    const Transform& o2w, const Transform& w2o, bool reverseOrientation,
    const int nTriangles, const int nVertices, const vector<Point3f>& vPos, vector<int>& faceVID, vector<shared_ptr<material>>& m) {
    TriangleMeshNew mesh(o2w, reverseOrientation, nTriangles, nVertices, vPos, faceVID,m);
    vector<shared_ptr<Primitive>> meshT;
    //meshT.resize(nTriangles);

    for (int i = 0; i < nTriangles; i++) {
        Point3f fv[3] = { mesh.vPos[faceVID[3 * i]],mesh.vPos[faceVID[3 * i + 1]] ,mesh.vPos[faceVID[3 * i + 2]] };
        meshT.push_back(make_shared<Trianglenew >(o2w, w2o,
            mesh.reverseOrientation, fv, mesh.mat_ptr[i]));
    }

    return meshT;
}

void ReadObj(const char* filename, vector<Point3f>& vPos, vector<int>& faceVID, double scale) {
    FILE* file;
    file = fopen(filename, "r");

    Float a, b, c;
    int e, f, g;
    char v;

    vPos.clear();
    faceVID.clear();

    double xmax = -infinity, xmin = infinity, ymax = -infinity, ymin = infinity, zmax = -infinity, zmin = infinity;

    while (!feof(file))
    {
        v = fgetc(file);

        if (v == 'v')
        {
            fscanf(file, "%f%f%f", &a, &b, &c);
            vPos.push_back(Point3f(a, b, c));

            if (a > xmax) xmax = a;
            if (a < xmin) xmin = a;
            if (b > ymax) ymax = b;
            if (b < ymin) ymin = b;
            if (c > zmax) zmax = c;
            if (c < zmin) zmin = c;
            //cout << a << " " << b << " " << c << endl;
            
        }
        else if (v == 'f')
        {
            fscanf(file, "%d%d%d", &e, &f, &g);
            faceVID.push_back(e-1);
            faceVID.push_back(f-1);
            faceVID.push_back(g-1);
        }
    }


    //cout << xmax << " " << xmin << " " << ymax << " " << ymin << " " << zmax << " " << zmin << endl;
    // 
    // 
    //cout << vPos.size() << " " << faceVID.size()<<endl;

    for (int i = 0; i < vPos.size(); i++)
        vPos[i] *= scale;
}

std::vector<std::shared_ptr<Primitive>> CreateTriangleMesh(
    const Transform& o2w, const Transform& w2o, bool reverseOrientation,
    const int nTriangles, const int nVertices, const vector<Point3f>& vPos, vector<int>& faceVID, shared_ptr<material> m) {
    //TriangleMeshNew mesh(o2w, reverseOrientation, nTriangles, nVertices, vPos, faceVID, m);
    vector<shared_ptr<Primitive>> meshT;
    //meshT.resize(nTriangles);

    for (int i = 0; i < nTriangles; i++) {
        Point3f fv[3] = { ConvertPTrans(o2w,vPos[faceVID[3 * i]]),ConvertPTrans(o2w,vPos[faceVID[3 * i + 1]]) ,ConvertPTrans(o2w,vPos[faceVID[3 * i + 2]]) };
        //cout << fv[0].x << endl;
        meshT.push_back(make_shared<Trianglenew >(o2w, w2o,
            reverseOrientation, fv, m));
    }

    return meshT;
}

#ifdef USE_OPENMESH
void ReadObjOpenMesh(const char* filename, vector<Point3f>& vPos, vector<int>& faceVID, double scale) {
    Mesh mesh;
    OpenMesh::IO::read_mesh(mesh, filename);

    for (int i = 0; i < mesh.n_faces(); i++) {
        Mesh::FaceVertexIter fv_it = mesh.fv_begin(mesh.face_handle(i));
        for (; fv_it.is_valid(); fv_it++) {
            faceVID.push_back((*fv_it).idx());
        }
    }
    for (int i = 0; i < mesh.n_vertices(); i++) {
        float* p = mesh.point(mesh.vertex_handle(i)).data();
        vPos.push_back(Point3f(p[0]*scale,p[1] * scale,p[2] * scale));
    }
    cout << vPos.size() << " " << faceVID.size() << endl;
}
#endif // USE_OPENMESH
