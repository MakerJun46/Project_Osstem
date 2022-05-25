#include "stdafx.h"

#include "MyTriMesh.h"
#include "igl/exact_geodesic.h"
#include "igl/readOBJ.h"
#include "../Utility/SimpleCompute.h"
#include "quadricfitting.h"
#include "../MeshDeform/PoissonDeform/PoissonDeformation.h"

Vector3D g_zeroVec(0, 0, 0);


MyTriMesh::MyTriMesh(void)
{
	init();
}

MyTriMesh::~MyTriMesh(void)
{
	m_mesh.clear();
}

void MyTriMesh::init()
{
	m_disType = BOUNDARY_MESH;
	m_show = true;
	m_drawMesh = 1;
	m_updateDistList = true;
	boundBoxEdgeLen = 1;
}

void MyTriMesh::SetMesh(TriMesh* mesh_)
{
	m_mesh = *mesh_;
	UpdateDisList();
};

void  MyTriMesh::DrawList()
{
	if (m_mesh.vertices_empty())
		return;

	if (m_updateDistList)
		m_updateDistList = false;
	else
		return;

	//create display list
	glNewList(m_drawMesh, GL_COMPILE);

	TriMesh::FaceIter f_it;
	TriMesh::FaceVertexIter fv_it;

	//move object to origin
	glTranslated(-m_centerOfAABB.m_x, -m_centerOfAABB.m_y, -m_centerOfAABB.m_z);

	//draw the model with GL_LINE_LOOP
	glLineWidth(1);


	int colorIndex;

	//怜삥齡듐
	if (m_disType == POINT_ONLY_MESH)
	{
		glDisable(GL_LIGHTING);
		glColor3f(1, 0, 0);
		glBegin(GL_POINTS);
		glPointSize(5);
		for (TriMesh::VertexIter vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); vi++)
		{
			break;
			TriMesh::Point point = m_mesh.point(vi);
			glVertex3d(point[0], point[1], point[2]);
		}
		glEnd();
		glPointSize(1);
		glEndList();
		return;
	}

	//삥齡窟움
	if (m_disType != SMOOTH_MESH)
	{
		glColor3f(0.2, 0.1, 0.2);
		glDisable(GL_LIGHTING);

		for (f_it = m_mesh.faces_begin(); f_it != m_mesh.faces_end(); f_it++)
		{
			glBegin(GL_LINE_LOOP);
			for (fv_it = m_mesh.fv_begin(f_it); fv_it != m_mesh.fv_end(f_it); fv_it++)
			{
				glVertex3d(m_mesh.point(fv_it)[0], m_mesh.point(fv_it)[1], m_mesh.point(fv_it)[2]);
			}
			glEnd();
		}
	}

	//삥齡밟뺄충
	{
		//get mormal
		m_mesh.request_face_normals();
		m_mesh.request_vertex_normals();
		m_mesh.update_normals();

		//offset z-depth
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1, 1);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//show the material visual color
		glDisable(GL_COLOR_MATERIAL);

		//draw the model with GL_POLYGON
		glEnable(GL_LIGHTING);
		for (f_it = m_mesh.faces_begin(); f_it != m_mesh.faces_end(); f_it++)
		{
			glBegin(GL_POLYGON);
			glNormal3d(m_mesh.normal(f_it)[0], m_mesh.normal(f_it)[1], m_mesh.normal(f_it)[2]);
			for (fv_it = m_mesh.fv_begin(f_it); fv_it != m_mesh.fv_end(f_it); fv_it++)
			{
				glNormal3d(m_mesh.normal(fv_it)[0], m_mesh.normal(fv_it)[1], m_mesh.normal(fv_it)[2]);
				//glNormal3d(m_vertexNormal[fv_it.handle().idx()].m_x,m_vertexNormal[fv_it.handle().idx()].m_y,m_vertexNormal[fv_it.handle().idx()].m_z);
				glVertex3d(m_mesh.point(fv_it)[0], m_mesh.point(fv_it)[1], m_mesh.point(fv_it)[2]);
			}
			glEnd();
		}
		glPopAttrib();
	}
	glEndList();
}

//뵙懃삥齡변鑒
void MyTriMesh::Draw()
{
	if (!m_show)
		return;
	DrawList();
	glCallList(m_drawMesh);
	glFlush();
}


void MyTriMesh::GetAABBofObject()
{
	if (m_mesh.vertices_empty())
		return;

	TriMesh::VertexIter vbegin = m_mesh.vertices_begin();
	TriMesh::Point pbegin = m_mesh.point(vbegin);

	double m_leastX = pbegin[0], m_mostX = pbegin[0];
	double m_leastY = pbegin[1], m_mostY = pbegin[1];
	double m_leastZ = pbegin[2], m_mostZ = pbegin[2];

	TriMesh::Point vertex;
	for (TriMesh::VertexIter vi = vbegin; vi != m_mesh.vertices_end(); vi++)
	{
		vertex = m_mesh.point(vi);

		//minx & maxx
		if (vertex[0] < m_leastX)
			m_leastX = vertex[0];
		else if (m_mostX < vertex[0])
			m_mostX = vertex[0];

		//miny & maxy
		if (vertex[1] < m_leastY)
			m_leastY = vertex[1];
		else if (m_mostY < vertex[1])
			m_mostY = vertex[1];

		//minz & maxz
		if (vertex[2] < m_leastZ)
			m_leastZ = vertex[2];
		else if (m_mostZ < vertex[2])
			m_mostZ = vertex[2];
	}

	m_centerOfAABB.m_x = (m_leastX + m_mostX) / 2.0;
	m_centerOfAABB.m_y = (m_leastY + m_mostY) / 2.0;
	m_centerOfAABB.m_z = (m_leastZ + m_mostZ) / 2.0;

	boundBoxEdgeLen = m_mostX - m_leastX;
	if (m_mostY - m_leastY > boundBoxEdgeLen)
		boundBoxEdgeLen = m_mostY - m_leastY;
	if (m_mostZ - m_leastZ > boundBoxEdgeLen)
		boundBoxEdgeLen = m_mostZ - m_leastZ;
}

double MyTriMesh::GetAverageEdgeLength()
{
	double averageEdgeLen = 0;
	int numEdges = 0;
	for (TriMesh::VertexIter v_it = m_mesh.vertices_begin(); v_it != m_mesh.vertices_end(); ++v_it)
	{
		TriMesh::Point point0 = m_mesh.point(v_it);
		for (TriMesh::VertexVertexIter vv_it = m_mesh.vv_iter(v_it); vv_it; ++vv_it)
		{
			TriMesh::Point point1 = m_mesh.point(vv_it);
			averageEdgeLen += (point0 - point1).length();
			numEdges++;
		}
	}
	averageEdgeLen /= 2 * numEdges;
	return averageEdgeLen;
}


void MyTriMesh::ReComputeVertexNormal()
{
	m_mesh.request_face_normals();
	m_mesh.request_vertex_normals();
	m_mesh.update_normals();

	m_vertexNormal.resize(m_mesh.n_vertices());

	for (TriMesh::VertexIter vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); vi++)
	{
		double sumArea = 0;
		vector<double> triAreas;
		vector<Vector3D> triFaceNorm;

		for (TriMesh::VertexFaceIter vf_i = m_mesh.vf_begin(vi); vf_i != m_mesh.vf_end(vi); vf_i++)
		{
			int i = 0;
			Point3D triPoint[3];
			triFaceNorm.push_back(Vector3D(m_mesh.normal(vf_i)[0], m_mesh.normal(vf_i)[1], m_mesh.normal(vf_i)[2]));

			for (TriMesh::FaceVertexIter fv_i = m_mesh.fv_begin(vf_i); fv_i != m_mesh.fv_end(vf_i); fv_i++)
			{
				TriMesh::Point point = m_mesh.point(fv_i);
				triPoint[i++] = Point3D(point[0], point[1], point[2]);
			}

			double area = SimpleCompute::GetTriangleArea(triPoint[0], triPoint[1], triPoint[2]);
			triAreas.push_back(area);
			sumArea += area;
		}

		Vector3D avgNorm;
		for (int i = 0; i < triAreas.size(); i++)
		{
			double rato = triAreas[i] / sumArea;
			avgNorm += triFaceNorm[i] * rato;
		}

		avgNorm.mf_normalize();
		m_vertexNormal[vi.handle().idx()] = avgNorm;
	}
}


void MyTriMesh::SaveMesh()
{
	if (m_mesh.vertices_empty())
		return;

	CString m_filePath;
	CFileDialog file(false);
	file.m_ofn.lpstrFilter = _T("mesh file(*.obj)\0*.obj\0Off File(*.off)\0*.off\0All File(*.*)\0*.*\0\0");
	file.m_ofn.lpstrTitle = _T("SAVE");
	if (file.DoModal() == IDOK)
	{
		m_filePath = file.GetPathName();
		OpenMesh::IO::write_mesh(m_mesh, m_filePath.GetBuffer(0));
	}
}

//==================================================OpenMesh=========================
#include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
 // -------------------- OpenMesh

 //--------------------- allquadrics

#include <fstream>
#include <map>
 //--------------------- allquadrics
#include <numeric>

#ifdef _DEBUG

#define debug_msg(msg) cout << msg <<endl
#define loop_confirm() debug_msg("right loop confirm")
#define loop_abnormal() debug_msg("Error : loop abnormal")

#else

#define debug_msg(msg)
#define loop_confirm()
#define loop_abnormal()

#endif

using namespace allquadrics;
using namespace std;

extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }
// ---------------------------------Calc Coordinate-------------------------------------------

#define SMALL_NUMMBER 0.000001

#define IS_ZERO(X) (X > -SMALL_NUMMBER && X < SMALL_NUMMBER)

#define IS_POSITIVE(X) (X > SMALL_NUMMBER)
#define IS_NEGATIVE(X) (X < -SMALL_NUMMBER)

vector<double> calc_abc(vec3 n, vec3 point, allquadrics::Quadric qfit)
{
	// a => anx^2 + bny^2 + cnz^2 + dnxny + enxnz + fnynz
	double a = (qfit.q[4] * n[0] * n[0]) + (qfit.q[7] * n[1] * n[1]) + (qfit.q[9] * n[2] * n[2])
		+ (qfit.q[5] * n[0] * n[1]) + (qfit.q[6] * n[0] * n[2]) + (qfit.q[8] * n[1] * n[2]);

	// b => 2axnx + 2byny + 2cznz + dnxy + dxny + enxz + exnz + fnyz + fynz + hnx + iny + jnz
	double b = (2 * qfit.q[4] * point[0] * n[0]) + (2 * qfit.q[7] * point[1] * n[1]) + (2 * qfit.q[9] * point[2] * n[2])
		+ (qfit.q[5] * n[0] * point[1]) + (qfit.q[5] * point[0] * n[1])
		+ (qfit.q[6] * n[0] * point[2]) + (qfit.q[6] * point[0] * n[2])
		+ (qfit.q[8] * n[1] * point[2]) + (qfit.q[8] * point[1] * n[2])
		+ (qfit.q[1] * n[0]) + (qfit.q[2] * n[1]) + (qfit.q[3] * n[2]);

	// c => ax^2 + by^2 + cz^2 + dxy + exz + fyz + hx + iy + jz + k
	double c = (qfit.q[4] * point[0] * point[0]) + (qfit.q[7] * point[1] * point[1]) + (qfit.q[9] * point[2] * point[2])
		+ (qfit.q[5] * point[0] * point[1]) + (qfit.q[6] * point[0] * point[2]) + (qfit.q[8] * point[1] * point[2])
		+ (qfit.q[1] * point[0]) + (qfit.q[2] * point[1]) + (qfit.q[3] * point[2]) + qfit.q[0];

	return { a, b, c };

}

vec3 partial_derivative(allquadrics::Quadric qfit, vec3 point)
{
	//ax^2 + by^2 + cz^2 + dxy + exz + fyz + hx + iy + jz + k = 0;

				 //  2 a x     +    d y    +    ez     +    h
	double x = (2 * qfit.q[4] * point[0]) + qfit.q[5] * point[1] + qfit.q[6] * point[2] + qfit.q[1];
	//  2 b y     +    d x    +    fz     +    i
	double y = (2 * qfit.q[7] * point[1]) + qfit.q[5] * point[0] + qfit.q[8] * point[2] + qfit.q[2];
	//  2 c z     +    e x    +    fy     +    j
	double z = (2 * qfit.q[9] * point[2]) + qfit.q[6] * point[0] + qfit.q[8] * point[1] + qfit.q[3];

	return vec3(x, y, z);	// (nx, ny, nz)
}

vector<double> solve_quadratic_equation(double a, double b, double c)
{
	vector<double> return_value;

	if (IS_ZERO(a))
	{
		if (IS_ZERO(b))
		{
			if (IS_ZERO(c))
			{
				for (int i = 0; i < 3; i++)
				{
					return_value.push_back(0);
				}
			}
		}
		else
		{
			return_value.push_back({ (-c / b) });
		}
	}
	else // 근의 공식
	{
		if (IS_POSITIVE((b * b) - (4 * a * c)))
		{
			double x = (-b + sqrt((b * b) - (4 * a * c))) / (2 * a);

			return_value.push_back(x);

			x = (-b - sqrt((b * b) - (4 * a * c))) / (2 * a);

			return_value.push_back(x);
		}
		else if (IS_ZERO((b * b) - (4 * a * c)))
		{
			return_value.push_back((-b / (2 * a)));
		}
	}

	return return_value;
}

double quadrics_equation(allquadrics::Quadric qfit, double x, double y, double z)
{
	return qfit.q[4] * x * x + qfit.q[7] * y * y + qfit.q[9] * z * z
		+ qfit.q[5] * x * y + qfit.q[6] * x * z + qfit.q[8] * y * z
		+ qfit.q[1] * x + qfit.q[2] * y + qfit.q[3] * z + qfit.q[0];
}

vec3 calc_target_point(allquadrics::Quadric qfit, vec3 point)
{
	double solution = quadrics_equation(qfit, point[0], point[1], point[2]);

	vec3 pd = partial_derivative(qfit, point);	// nx, ny, nz
	vector<double> abc = calc_abc(pd, point, qfit); // a, b, c

	vector<double> solution_quadratic = solve_quadratic_equation(abc[0], abc[1], abc[2]);  // t

	vector<double> positive_num;
	vector<double> negative_num;

	for (int i = 0; i < solution_quadratic.size(); i++)
	{
		if (IS_POSITIVE(solution_quadratic[i]))
		{
			positive_num.push_back(solution_quadratic[i]);
		}
		else
		{
			negative_num.push_back(solution_quadratic[i]);
		}
	}

	sort(positive_num.begin(), positive_num.end());	// 0번째가 positive 중 가장 작은 것
	sort(negative_num.rbegin(), negative_num.rend()); // 0번째가 negative 중 가장 큰 것

	//p(t) = (x y z) + (nx ny nz) * t
	if (IS_POSITIVE(solution) && !negative_num.empty())
	{
		double x = point[0] + pd[0] * negative_num[0];// x + nx * t
		double y = point[1] + pd[1] * negative_num[0];// y + ny * t
		double z = point[2] + pd[2] * negative_num[0];// z + nz * t

		vec3 v(x, y, z);

		return v;
	}
	else if (IS_NEGATIVE(solution) && !positive_num.empty())
	{
		double x = point[0] + pd[0] * positive_num[0];// x + nx * t
		double y = point[1] + pd[1] * positive_num[0];// y + ny * t
		double z = point[2] + pd[2] * positive_num[0];// z + nz * t

		vec3 v(x, y, z);

		return v;
	}

	return point;
}
//-------------------------------------------------------------------------------------------------


#define M_PI 3.14159

struct UnitHole
{
	int hole_number;
	bool is_filled;
	vector<pair<TriMesh::VertexHandle, TriMesh::Normal>> original_vertices_normals;
	vector<pair<TriMesh::VertexHandle, TriMesh::Point>> hole_original_vertices;
	vector<TriMesh::FaceHandle> original_Faces;
	vector<TriMesh::VertexHandle> new_Vertices;
	vector<TriMesh::FaceHandle> new_Faces;
	vector<pair<int, double>> vertexIndex_angles;
};

void save(TriMesh mesh, string name, bool resultSave = false)
{
	try
	{
		if (resultSave)
		{
			if (!OpenMesh::IO::write_mesh(mesh, name + ".obj"))
			{
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			}
		}
		else
		{
			if (!OpenMesh::IO::write_mesh(mesh, name + ".obj"))
			{
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			}
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}

bool cmpSecond(pair<int, double> a, pair<int, double> b)
{
	return a.second < b.second;
}

// calc edges average length
double edge_length_average(TriMesh& mesh)
{
	int count = 0;
	double edge_length_average = 0;

	for (auto it = mesh.edges_begin(); it != mesh.edges_end(); it++)
	{
		edge_length_average += mesh.calc_edge_length(*it);
		count++;
	}

	return edge_length_average /= count;
}

// set Original Vertices, Face handles
void set_isOriginal_Vertex_array(TriMesh& mesh, map<int, bool>& isOriginal_Vertex)
{
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		isOriginal_Vertex[it->idx()] = true;
	}
}
void set_isOriginal_Face_array(TriMesh& mesh, map<int, bool>& isOriginal_Face)
{
	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		isOriginal_Face[it->idx()] = true;
	}
}

// find boundary vertices and sort loop
void boundary_vertex_Sort(TriMesh mesh, vector<TriMesh::VertexHandle> vhandle_boundary, vector<pair<TriMesh::VertexHandle, TriMesh::Point>>& holes, vector<TriMesh::VertexHandle>& second_hole)
{
	holes.clear();

	if (vhandle_boundary.empty())
		return;

	holes.push_back(make_pair(vhandle_boundary[0], mesh.point(vhandle_boundary[0])));
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	TriMesh::HalfedgeHandle connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
	TriMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfEdge);
	vhandle_boundary.erase(vhandle_boundary.begin());

	while (!vhandle_boundary.empty())
	{
		if (vhandle_boundary.size() == 1)
		{
			vertex_tmp = vhandle_boundary[0];
			if (vertex_tmp.idx() == next_v.idx())
			{
				holes.push_back(make_pair(vertex_tmp, mesh.point(vertex_tmp)));

				connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
				next_v = mesh.to_vertex_handle(connected_halfEdge);


				if (next_v.idx() == vertex_startPoint.idx())
				{
					loop_confirm();
				}
				else
				{
					loop_abnormal();
					break;
				}
				return;
			}
			else
			{
				auto p = mesh.point(vertex_tmp);

				loop_abnormal();

				break;
			}
		}

		auto target_it = find(vhandle_boundary.begin(), vhandle_boundary.end(), next_v);

		if (target_it == vhandle_boundary.end() && next_v == vertex_startPoint)
		{
			//cout << "구멍이 2개 이상임, 남은 vertex들은 다른 구멍의 vertex" << endl;
			second_hole = vhandle_boundary; // 남은 다른 구멍의 vertex들을 따로 저장
			vhandle_boundary.clear();
		}
		else
		{
			holes.push_back(make_pair(*target_it, mesh.point(*target_it)));
			vertex_tmp = *target_it;

			connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
			next_v = mesh.to_vertex_handle(connected_halfEdge);

			vhandle_boundary.erase(target_it);
		}
	}
}
void boundary_Update(TriMesh& mesh, vector<pair<TriMesh::VertexHandle, TriMesh::Point>>& holes, vector<TriMesh::VertexHandle>& second_hole)
{
	vector<TriMesh::VertexHandle> vhandle_boundary;

	if (holes.size() <= 2)	// 이 구멍의 vertex가 2개 이하인 경우 => 다채운 경우이므로 다음 구멍으로 넘어감
	{
		for (int i = 0; i < second_hole.size(); i++)
		{
			if (mesh.is_valid_handle(second_hole[i]) && mesh.is_boundary(second_hole[i]))
			{
				vhandle_boundary.push_back(second_hole[i]);
			}
		}

		second_hole.clear(); // 다 채워주었으면 clear
	}
	else // 2개 이상인 경우는 아직 이을 vertex가 남아있으므로 해당 구멍에서 검사
	{
		for (int i = 0; i < holes.size(); i++)
		{
			if (mesh.is_valid_handle(holes[i].first) && mesh.is_boundary(holes[i].first))
			{
				vhandle_boundary.push_back(holes[i].first);
			}
		}
	}

	boundary_vertex_Sort(mesh, vhandle_boundary, holes, second_hole);
}

double calc_sector_angle_UseQaudrics(TriMesh& mesh, TriMesh::HalfedgeHandle _in_heh, vec3 normal)
{
	TriMesh::Normal v0, v1;
	TriMesh::Normal f_n(normal[0], normal[1], normal[2]);
	mesh.calc_sector_vectors(_in_heh, v0, v1);

	TriMesh::Scalar a = norm(v0);
	TriMesh::Scalar b = norm(v1);

	TriMesh::Scalar denom = norm(v0) * norm(v1);
	if (denom == TriMesh::Scalar(0))
	{
		return 0;
	}
	TriMesh::Scalar cos_a = dot(v0, v1) / denom;
	if (mesh.is_boundary(_in_heh))
	{
		TriMesh::Scalar sign_a = dot(cross(v0, v1), f_n);
		return OpenMesh::angle(cos_a, sign_a);
	}
	else
	{
		if (cos_a < -1)
			cos_a = -1;
		if (cos_a > 1)
			cos_a = 1;

		return acos(cos_a);
	}
}

// calc vertices angles and sort
double calc_angle(TriMesh& mesh, TriMesh::VertexHandle vh, bool useQuadrics, allquadrics::Quadric qfit)
{
	int idx = vh.idx();

	double deg;

	////if (useQuadrics)
	////{
	//MyMesh::Point p = mesh.point(vh);

	//vec3 point(p[0], p[1], p[2]);

	//vec3 normal = partial_derivative(qfit, point);

	//MyMesh::HalfedgeHandle prev_halfedge = mesh.prev_halfedge_handle(mesh.halfedge_handle(vh));

	//deg = calc_sector_angle_UseQaudrics(mesh, prev_halfedge, normal);
	//}
	//else
	//{
	TriMesh::HalfedgeHandle prev_halfedge = mesh.prev_halfedge_handle(mesh.halfedge_handle(vh));

	deg = mesh.calc_sector_angle(prev_halfedge) * 180 / M_PI;
	//}

	if (deg < 0)
		deg += 360;

	return deg;
}

void angle_Update(TriMesh& mesh, vector<pair<TriMesh::VertexHandle, TriMesh::Point>>& holes, vector<pair<int, double>>& index_angles, bool useQaudrics, allquadrics::Quadric qfit)
{
	index_angles.clear();
	mesh.request_face_normals();

	for (int i = 0; i < holes.size(); i++)
	{
		index_angles.push_back(make_pair(i, calc_angle(mesh, holes[i].first, useQaudrics, qfit)));
	}

	sort(index_angles.begin(), index_angles.end(), cmpSecond);
}

// merge close vertices
bool isClose_OtherVertex_Exist(vector<pair<TriMesh::VertexHandle, TriMesh::Point>>& holes, int a, int b)
{
	TriMesh::Point va1; // a 이전 vertex Point
	TriMesh::Point va2; // a 다음 vertex Point
	TriMesh::Point vb1; // b 이전 vertex Point
	TriMesh::Point vb2; // b 다음 vertex Point

	TriMesh::Point va = holes[a].second;
	TriMesh::Point vb = holes[b].second;

	double interval = sqrt(pow((va[0] - vb[0]), 2) + pow((va[1] - vb[1]), 2) + pow((va[2] - vb[2]), 2));

	va1 = a == 0 ? holes[holes.size() - 1].second : holes[a - 1].second;
	va2 = a == holes.size() - 1 ? holes[0].second : holes[a + 1].second;
	vb1 = b == 0 ? holes[holes.size() - 1].second : holes[b - 1].second;
	vb2 = b == holes.size() - 1 ? holes[0].second : holes[b + 1].second;

	vector<double> cmp_intervals;

	cmp_intervals.push_back(sqrt(pow((va[0] - vb1[0]), 2) + pow((va[1] - vb1[1]), 2) + pow((va[2] - vb1[2]), 2)));
	cmp_intervals.push_back(sqrt(pow((va[0] - vb2[0]), 2) + pow((va[1] - vb2[1]), 2) + pow((va[2] - vb2[2]), 2)));
	cmp_intervals.push_back(sqrt(pow((vb[0] - va1[0]), 2) + pow((vb[1] - va1[1]), 2) + pow((vb[2] - va1[2]), 2)));
	cmp_intervals.push_back(sqrt(pow((vb[0] - va2[0]), 2) + pow((vb[1] - va2[1]), 2) + pow((vb[2] - va2[2]), 2)));

	for (double i : cmp_intervals)
	{
		if (i <= interval)
		{
			return true;
		}
	}

	return false;
}
void fill_triangles_afterMerge(TriMesh& mesh, TriMesh::VertexHandle currentVh)
{
	TriMesh::HalfedgeHandle prev_halfedge = mesh.halfedge_handle(currentVh);

	vector<TriMesh::VertexHandle> TrianglePoints;

	int count = 0;

	while (count < 2)
	{
		prev_halfedge = mesh.prev_halfedge_handle(prev_halfedge);
		TriMesh::VertexHandle targetVh = mesh.from_vertex_handle(prev_halfedge);

		if (targetVh == currentVh)
		{
			TrianglePoints.push_back(mesh.to_vertex_handle(prev_halfedge));
			TrianglePoints.push_back(mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(currentVh))));
			TrianglePoints.push_back(currentVh);

			mesh.add_face(TrianglePoints);

			TrianglePoints.clear();

			prev_halfedge = mesh.halfedge_handle(currentVh);

			count++;
		}
	}
}
void merge_close_vertex(TriMesh& mesh, vector<pair<TriMesh::VertexHandle, TriMesh::Point>>& holes, map<int, bool> isOriginal_Vertex, double edge_length_average, int k, bool useQuadrics, allquadrics::Quadric qfit)
{
	for (int i = 0; i < holes.size(); i++)
	{
		vec3 p1(holes[i].second[0], holes[i].second[1], holes[i].second[2]);

		for (int j = 0; j < holes.size(); j++)
		{
			if (holes.size() > 5)
			{
				if (((i - j + holes.size()) % holes.size() <= 2) || ((j - i + holes.size()) % holes.size() <= 2))
					continue;
			}
			else
			{
				if (((i - j + holes.size()) % holes.size() <= 1) || ((j - i + holes.size()) % holes.size() <= 1))
					continue;
			}

			vec3 p2(holes[j].second[0], holes[j].second[1], holes[j].second[2]);

			double interval = sqrt(pow((p1[0] - p2[0]), 2) + pow((p1[1] - p2[1]), 2) + pow((p1[2] - p2[2]), 2));

			TriMesh::HalfedgeHandle hh1 = mesh.find_halfedge(holes[i].first, holes[j].first);
			TriMesh::HalfedgeHandle hh2 = mesh.find_halfedge(holes[j].first, holes[i].first);

			if (interval < edge_length_average && !isOriginal_Vertex[holes[j].first.idx()] && !isOriginal_Vertex[holes[i].first.idx()]
				&& hh1 == TriMesh::InvalidHalfedgeHandle && hh2 == TriMesh::InvalidHalfedgeHandle && !isClose_OtherVertex_Exist(holes, i, j))
			{
				vec3 center_coord((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2);

				if (useQuadrics)
					center_coord = calc_target_point(qfit, center_coord);

				TriMesh::Point p_new(center_coord[0], center_coord[1], center_coord[2]);

				mesh.set_point(holes[j].first, p_new);

				vector<vector<TriMesh::VertexHandle>> triangle_Points;

				for (auto it = mesh.vf_begin(holes[i].first); it != mesh.vf_end(holes[i].first); it++)	// merge vertices
				{
					TriMesh::FaceHandle fh = it;
					vector<TriMesh::VertexHandle> tmp;
					for (auto it_f = mesh.fv_begin(fh); it_f != mesh.fv_end(fh); it_f++)
					{
						if (*it_f == holes[i].first)
							tmp.push_back(holes[j].first);
						else
							tmp.push_back(it_f);
					}
					triangle_Points.push_back(tmp);
					tmp.clear();
				}

				mesh.request_face_status();
				mesh.request_edge_status();
				mesh.request_vertex_status();

				auto it = holes.begin() + i;
				TriMesh::VertexHandle target = holes[i].first;
				TriMesh::VertexHandle merge_vh = holes[j].first;

				holes.erase(it);

				mesh.delete_vertex(target);

				for (auto vh : triangle_Points)
				{
					vector<TriMesh::VertexHandle> tmp;
					for (auto v : vh)
					{
						tmp.push_back(v);
					}

					mesh.add_face(tmp);
					tmp.clear();
				}

				fill_triangles_afterMerge(mesh, merge_vh);

				mesh.garbage_collection();

				mesh.request_face_status();
				mesh.request_edge_status();
				mesh.request_vertex_status();

				//save(mesh, "after_merge_" + to_string(k));

				return;
			}
		}
	}
}

// fill hole
void fill_hole(TriMesh& mesh_findBoundary, vector<pair<TriMesh::VertexHandle, TriMesh::Point>> holes, vector<pair<int, double>> index_angles,
	map<int, bool> isOriginal_Vertex, double edge_length_average, allquadrics::Quadric qfit, bool useQuadrics)
{
	vector<vector<TriMesh::VertexHandle>> triangle_Points;
	vector<TriMesh::VertexHandle> second_hole;

	int k = 0;

	while (!holes.empty())
	{
		string msg = "동작 중 : " + to_string(k) + "번째 반복 실행";
		debug_msg(msg);

		//for (auto i : index_angles)
		//{
		//	cout << i.first << "번의 사이각 : " << i.second << endl;
		//}

		if (k % 100 == 0)
			save(mesh_findBoundary, to_string(k));

		if (holes.size() < 3)
		{
			break;
		}



		int index = index_angles[0].first;	// 사이각이 가장 작은 개체(내림차순 정렬의 0번째)
		vector<TriMesh::VertexHandle> tmp;
		vec3 index_vec(holes[index].second[0], holes[index].second[1], holes[index].second[2]);
		bool newVertexCreated;


		if (index_angles[0].second <= 75.1) // 75보다 각이 작은 경우 -> 점 0개 추가, 점 1개 제외, 삼각형 1개 추가
		{
			if (index == 0)
			{
				tmp.push_back(holes[holes.size() - 1].first);
			}
			else
			{
				tmp.push_back(holes[index - 1].first);
			}

			tmp.push_back(holes[index].first);

			if (index == holes.size() - 1)
			{
				tmp.push_back(holes[0].first);
			}
			else
			{
				tmp.push_back(holes[index + 1].first);
			}

			triangle_Points.push_back(tmp); // 삼각형 추가

			holes.erase(find(holes.begin(), holes.end(), holes[index]));

			newVertexCreated = false;
		}

		else if (75.1 < index_angles[0].second && index_angles[0].second <= 135.1) // 75 < 사이각 < 135 인 경우 -> 점 1개 추가, 점 1개 제외, 삼각형 2개 추가
		{
			vec3 a;
			if (index == 0)
			{
				a[0] = holes[holes.size() - 1].second[0];
				a[1] = holes[holes.size() - 1].second[1];
				a[2] = holes[holes.size() - 1].second[2];

				tmp.push_back(holes[holes.size() - 1].first);
			}
			else
			{
				a[0] = holes[index - 1].second[0];
				a[1] = holes[index - 1].second[1];
				a[2] = holes[index - 1].second[2];

				tmp.push_back(holes[index - 1].first);
			}

			vec3 b;
			if (index == holes.size() - 1)
			{
				b[0] = holes[0].second[0];
				b[1] = holes[0].second[1];
				b[2] = holes[0].second[2];
			}
			else
			{
				b[0] = holes[index + 1].second[0];
				b[1] = holes[index + 1].second[1];
				b[2] = holes[index + 1].second[2];
			}
			vec3 v1 = a - index_vec;
			vec3 v2 = b - index_vec;
			vec3 new_vec = v1 + v2; // 새로운 점 좌표 (a 벡터 + b 벡터)

			new_vec.normalize();

			new_vec = new_vec;// *((v1.length() + v2.length()) / 2);

			new_vec *= max(edge_length_average, ((v1.length() + v2.length()) / 2));

			vec3 new_coord = index_vec + new_vec;

			if (useQuadrics && holes.size() > 10)
				new_coord = calc_target_point(qfit, new_coord);	// calc quadrics

			TriMesh::Point new_point(new_coord[0], new_coord[1], new_coord[2]);

			TriMesh::VertexHandle vh_tmp = mesh_findBoundary.add_vertex(new_point); // 새로운 점 생성

			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp);
			triangle_Points.push_back(tmp); // 삼각형 1 추가

			tmp.clear();
			tmp.push_back(vh_tmp);
			tmp.push_back(holes[index].first);

			if (index == holes.size() - 1)
				tmp.push_back(holes[0].first);
			else
				tmp.push_back(holes[index + 1].first);

			triangle_Points.push_back(tmp); // 삼각형 2 추가

			holes[index] = make_pair(vh_tmp, new_point);// 삼각형 추가했으니 holes와 angle의 정보 새 점으로 바꾸기

			newVertexCreated = true;
		}

		else // 135 < 사이각 인 경우 -> 점 2개 추가, 점 1개 제외, 삼각형 3개 추가
		{
			vec3 a;
			if (index == 0)
			{
				a[0] = holes[holes.size() - 1].second[0];
				a[1] = holes[holes.size() - 1].second[1];
				a[2] = holes[holes.size() - 1].second[2];

				tmp.push_back(holes[holes.size() - 1].first);
			}
			else
			{
				a[0] = holes[index - 1].second[0];
				a[1] = holes[index - 1].second[1];
				a[2] = holes[index - 1].second[2];

				tmp.push_back(holes[index - 1].first);
			}
			vec3 b;
			if (index == holes.size() - 1)
			{
				b[0] = holes[0].second[0];
				b[1] = holes[0].second[1];
				b[2] = holes[0].second[2];
			}
			else
			{
				b[0] = holes[index + 1].second[0];
				b[1] = holes[index + 1].second[1];
				b[2] = holes[index + 1].second[2];
			}

			vec3 v1 = a - index_vec;
			vec3 v2 = b - index_vec;
			vec3 center_coord = v1 + v2;

			center_coord.normalize();

			center_coord = center_coord * ((v1.length() + v2.length()) / 2); //center coord = v1과 v2의 중점

			vec3 a1 = (v1 + center_coord * 2).normalize();// *(v1.length() + center_coord.length() * 2) / 3;
			vec3 b1 = (center_coord * 2 + v2).normalize();// *(center_coord.length() * 2 + v2.length()) / 3;
			a1 *= max(edge_length_average, (v1.length() + center_coord.length() * 2) / 3);
			b1 *= max(edge_length_average, (center_coord.length() * 2 + v2.length()) / 3);

			a1 += index_vec;
			b1 += index_vec;

			if (useQuadrics && holes.size() > 10)
			{
				a1 = calc_target_point(qfit, a1);	// calc quadrics
				b1 = calc_target_point(qfit, b1);
			}

			TriMesh::Point new_Point_a(a1[0], a1[1], a1[2]);
			TriMesh::Point new_Point_b(b1[0], b1[1], b1[2]);

			TriMesh::VertexHandle vh_tmp_a = mesh_findBoundary.add_vertex(new_Point_a);
			TriMesh::VertexHandle vh_tmp_b = mesh_findBoundary.add_vertex(new_Point_b);

			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp_a);
			triangle_Points.push_back(tmp);	// 삼각형 1 추가

			tmp.clear();
			tmp.push_back(vh_tmp_a);
			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp_b);
			triangle_Points.push_back(tmp); // 삼각형 2 추가

			tmp.clear();
			tmp.push_back(vh_tmp_b);
			tmp.push_back(holes[index].first);

			if (index == holes.size() - 1)
			{
				tmp.push_back(holes[0].first);
			}
			else
			{
				tmp.push_back(holes[index + 1].first);
			}

			triangle_Points.push_back(tmp);	// 삼각형 3 추가

			holes.insert(holes.begin() + index, make_pair(vh_tmp_b, new_Point_b)); // 새로운 point2 추가
			holes.insert(holes.begin() + index, make_pair(vh_tmp_a, new_Point_a));	// 새로운 point1 추가
			holes.erase(holes.begin() + index + 2); // point 삭제

			newVertexCreated = true;
		}

		//--------------------------------------------------- vertex들을 이어서 면으로 만들기
		for (int i = 0; i < triangle_Points.size(); i++)
		{
			vector<TriMesh::VertexHandle> face_h;

			for (int j = 0; j < triangle_Points[i].size(); j++)
			{
				face_h.push_back(triangle_Points[i][j]);
			}

			mesh_findBoundary.add_face(face_h);
		}
		//--------------------------------------------------- vertex들을 이어서 면으로 만들기


		triangle_Points.clear();
		tmp.clear();

		//save(mesh_findBoundary, to_string(k));
		if (newVertexCreated)
			merge_close_vertex(mesh_findBoundary, holes, isOriginal_Vertex, edge_length_average, k, useQuadrics, qfit);

		boundary_Update(mesh_findBoundary, holes, second_hole);
		angle_Update(mesh_findBoundary, holes, index_angles, useQuadrics, qfit);

		k++;
	}
}

// find new vertices, faces after fill hole process
bool is_other_hole_vertex(TriMesh::VertexHandle vh, int targethole_index, vector<UnitHole> unit_holes)
{
	for (int i = 0; i < unit_holes.size(); i++)
	{
		if (i == targethole_index)
			continue;

		for (auto v : unit_holes[i].new_Vertices)
		{
			if (vh == v)
				return true;
		}
	}

	return false;
}
bool is_other_hole_face(TriMesh::FaceHandle fh, int targethole_index, vector<UnitHole> unit_holes)
{
	for (int i = 0; i < unit_holes.size(); i++)
	{
		if (i == targethole_index)
			continue;

		for (auto f : unit_holes[i].new_Faces)
		{
			if (fh == f)
				return true;
		}
	}

	return false;
}
vector<TriMesh::VertexHandle> find_newVertices(TriMesh& mesh, UnitHole hole, vector<UnitHole> unit_holes, map<int, bool> isOriginal_Vertex)
{
	vector<TriMesh::VertexHandle> newVertexHandles;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		if (isOriginal_Vertex[it->idx()] || is_other_hole_vertex(it, hole.hole_number, unit_holes))
			continue;
		else
			newVertexHandles.push_back(*it);
	}

	return newVertexHandles;
}
vector<TriMesh::FaceHandle> find_newFaces(TriMesh& mesh, UnitHole hole, vector<UnitHole> unit_holes, map<int, bool> isOriginal_Face)
{
	vector<TriMesh::FaceHandle> newFaceHandles;

	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		if (isOriginal_Face[it->idx()] || is_other_hole_face(it, hole.hole_number, unit_holes))
			continue;
		else
			newFaceHandles.push_back(*it);
	}

	return newFaceHandles;
}

// find all holes
vector<TriMesh::VertexHandle> Original_boundary_vertex_Sort(TriMesh& mesh, vector<TriMesh::VertexHandle>& vhandle_boundary, vector<TriMesh::FaceHandle>& fhandle_boundary)
{
	vector<TriMesh::VertexHandle> hole_vertices;

	hole_vertices.push_back(vhandle_boundary[0]);
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	TriMesh::HalfedgeHandle connected_halfedge = mesh.halfedge_handle(vertex_tmp);
	TriMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfedge);
	vhandle_boundary.erase(vhandle_boundary.begin());

	while (true)
	{
		auto target_it = find(vhandle_boundary.begin(), vhandle_boundary.end(), next_v);

		fhandle_boundary.push_back(mesh.face_handle(mesh.opposite_halfedge_handle(connected_halfedge)));

		if (target_it == vhandle_boundary.end() && next_v == vertex_startPoint)
		{
			//cout << "구멍이 2개 이상임, 남은 vertex들은 다른 구멍의 vertex" << endl;
			return hole_vertices;
		}
		else
		{
			hole_vertices.push_back(*target_it);
			vertex_tmp = *target_it;


			connected_halfedge = mesh.halfedge_handle(vertex_tmp);
			next_v = mesh.to_vertex_handle(connected_halfedge);

			vhandle_boundary.erase(target_it);
		}

	}
}

void statusSettings(TriMesh& mesh)
{
	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();
	mesh.request_halfedge_status();

	mesh.update_normals();
}

void eraseFace_hasTwoHalfEdge_boundary(TriMesh& mesh)
{
	statusSettings(mesh);

	vector<TriMesh::FaceHandle> fh_tmp;

	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		int cnt = 0;

		for (auto it_he = mesh.fh_begin(it); it_he != mesh.fh_end(it); it_he++)
		{
			if (mesh.is_boundary(mesh.opposite_halfedge_handle(it_he)))
			{
				cnt++;
			}
		}

		if (cnt > 1)
		{
			cout << "EraseFace" << endl;
			fh_tmp.push_back(it);
		}
	}

	for (int i = 0; i < fh_tmp.size(); i++)
	{
		mesh.delete_face(fh_tmp[i]);
	}

	mesh.garbage_collection();

	save(mesh, "eraseTwoHalfedge_test");
}



vector<UnitHole> find_All_UnitHoles(TriMesh& mesh, bool useQuadrics, allquadrics::Quadric qfit)
{
	vector<UnitHole> vec_tmp;
	vector<TriMesh::VertexHandle> vhandle_boundary;
	vector<TriMesh::FaceHandle> fhandle_boundary;
	vector<TriMesh::VertexHandle> vhandle_tmp;
	int hole_number = 0;


	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		if (mesh.is_boundary(it))
			vhandle_boundary.push_back(*it);
	}

	while (!vhandle_boundary.empty())
	{
		vhandle_tmp = Original_boundary_vertex_Sort(mesh, vhandle_boundary, fhandle_boundary);
		UnitHole uh_tmp;

		for (auto vh : vhandle_tmp)
		{
			TriMesh::Normal vertex_normal = mesh.normal(vh);

			uh_tmp.original_vertices_normals.push_back(make_pair(vh, vertex_normal));
			uh_tmp.hole_original_vertices.push_back(make_pair(vh, mesh.point(vh)));
		}

		uh_tmp.original_Faces = fhandle_boundary;

		angle_Update(mesh, uh_tmp.hole_original_vertices, uh_tmp.vertexIndex_angles, useQuadrics, qfit);

		uh_tmp.hole_number = hole_number;
		uh_tmp.is_filled = false;

		vec_tmp.push_back(uh_tmp);

		cout << "찾은 구멍 index : " << hole_number << endl;

		hole_number++;
	}

	return vec_tmp;
}

void find_unvalid_Vertices(TriMesh& mesh, UnitHole& hole)
{
	mesh.request_vertex_status();
	mesh.request_face_status();

	for (int i = 0; i < hole.new_Vertices.size(); i++)
	{
		if (!mesh.is_valid_handle(hole.new_Vertices[i]) || mesh.is_boundary(hole.new_Vertices[i]))
		{
			mesh.delete_vertex(hole.new_Vertices[i]);
			hole.new_Vertices.erase(hole.new_Vertices.begin() + i);
			i--;

			mesh.garbage_collection();
		}
	}

}

vec3 cross_product_func(vec3 v1, vec3 v2)
{
	return vec3
	(
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]
	);
}

double dot_product_func(vec3 v1, vec3 v2)
{
	double sum = 0.0;

	sum += v1[0] * v2[0];
	sum += v1[1] * v2[1];
	sum += v1[2] * v2[2];

	return sum;
}

mat3 rotateAlign(vec3 v1, vec3 v2)
{
	vec3 axis = cross_product_func(v1, v2);

	const float cosA = dot_product_func(v1, v2);
	const float k = 1.0f / (1.0f + cosA);

	vec3 tmp((axis[0] * axis[0] * k) + cosA, (axis[1] * axis[0] * k) - axis[2], (axis[2] * axis[0] * k) + axis[1]);
	vec3 tmp1((axis[0] * axis[1] * k) + axis[2], (axis[1] * axis[1] * k) + cosA, (axis[2] * axis[1] * k) - axis[0]);
	vec3 tmp2((axis[0] * axis[2] * k) - axis[1], (axis[1] * axis[2] * k) + axis[0], (axis[2] * axis[2] * k) + cosA);

	/*mat3 result((axis[0] * axis[0] * k) + cosA,
		(axis[1] * axis[0] * k) - axis[2],
		(axis[2] * axis[0] * k) + axis[1],
		(axis[0] * axis[1] * k) + axis[2],
		(axis[1] * axis[1] * k) + cosA,
		(axis[2] * axis[1] * k) - axis[0],
		(axis[0] * axis[2] * k) - axis[1],
		(axis[1] * axis[2] * k) + axis[0],
		(axis[2] * axis[2] * k) + cosA
	);*/



	mat3 result(tmp, tmp1, tmp2);

	return result;
}

vec3 get_G(TriMesh::Point v1, TriMesh::Point v2, TriMesh::Point v3)
{
	return vec3((v1[0] + v2[0] + v3[0]) / 3,
		(v1[1] + v2[1] + v3[1] / 3),
		(v1[2] + v2[2] + v3[2] / 3));
}

void Deformation(TriMesh& mesh, UnitHole& hole, vector<bool>& isBoundaryVertex, VectorXd divMatri_x, VectorXd divMatri_y, VectorXd divMatri_z)
{
	MeshLaplaceSolver myLPLsolver;

	myLPLsolver.SetDesMesh(mesh);
	myLPLsolver.SetControlVertex(isBoundaryVertex);
	myLPLsolver.ComputeLalacianMatrixA();

	myLPLsolver.SetRightHandB(divMatri_x);
	VectorXd x = myLPLsolver.LplacianSolve();

	//solve y
	myLPLsolver.SetRightHandB(divMatri_y);
	VectorXd y = myLPLsolver.LplacianSolve();

	//solve z
	myLPLsolver.SetRightHandB(divMatri_z);
	VectorXd z = myLPLsolver.LplacianSolve();

	//update mesh
	TriMesh::VertexHandle vh;

	for (int i = 0; i < hole.new_Vertices.size(); i++)
	{
		vh = hole.new_Vertices[i];

		mesh.point(vh)[0] = x[vh.idx()];
		mesh.point(vh)[1] = y[vh.idx()];
		mesh.point(vh)[2] = z[vh.idx()];
	}


	//for (int i = 0; i < mesh.n_vertices(); i++)
	//{
	//	vh = mesh.vertex_handle(i);

	//	mesh.point(vh)[0] = x[i];
	//	mesh.point(vh)[1] = y[i];
	//	mesh.point(vh)[2] = z[i];
	//}
}

void Calc_geodesicDistance(TriMesh& mesh, UnitHole& hole, vector<vector<double>>& geodesicDistance)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::VectorXi VT, FT;
	Eigen::VectorXd d;

	vector<vector<double>> v_tmp;
	vector<vector<int>> f_tmp;
	vector<int> vt_tmp;
	vector<int> ft_tmp;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		TriMesh::Point p = mesh.point(it);

		vector<double> tmp;

		tmp.push_back(p[0]);
		tmp.push_back(p[1]);
		tmp.push_back(p[2]);

		v_tmp.push_back(tmp);
	}

	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		TriMesh::FaceHandle f = it;
		vector<int> tmp;

		for (auto f_it = mesh.fv_begin(f); f_it != mesh.fv_end(f); f_it++)
		{
			tmp.push_back(f_it->idx());
		}

		f_tmp.push_back(tmp);
	}

	for (int i = 0; i < hole.hole_original_vertices.size(); i++)
	{
		vt_tmp.push_back(hole.hole_original_vertices[i].first.idx());
	}

	for (int i = 0; i < hole.original_Faces.size(); i++)
	{
		ft_tmp.push_back(hole.original_Faces[i].idx());
	}

	igl::list_to_matrix(v_tmp, V);	// V matrix에 점의 좌표 입력 (행 : 점, 열 : 좌표값(double 3개))
	igl::list_to_matrix(f_tmp, F);  // F matrix에 face에 포함된 vertex의 idx 입력
	//igl::list_to_matrix(vt_tmp, VT);	// VT 에 boundary Verteices index 입력
	igl::list_to_matrix(ft_tmp, FT);	// FT 에 boundary Faces index 입력

	// new Vertices, new Faces를 하나씩 넣고 돌려서 distance 각각 저장

	for (int i = 0; i < hole.new_Faces.size(); i++)
	{
		TriMesh::FaceHandle fh = hole.new_Faces[i];

		Eigen::VectorXi VS, FS;
		VS.resize(0);
		FS.resize(1);

		auto it = mesh.fv_begin(fh);

		//VS[0] = it->idx();
		//it++;
		//VS[1] = it->idx();
		//it++;
		//VS[2] = it->idx();

		FS[0] = fh.idx();

		igl::exact_geodesic(V, F, VS, FS, VT, FT, d);


		vector<double> tmp;

		it = mesh.fv_begin(fh);

		for (int k = 0; k < d.size(); k++)
		{
			tmp.push_back(d[k]);
		}

		geodesicDistance[fh.idx()] = tmp; // geodesicDistance (face index) 마다 해당 Face의 모든 d 값 저장
	}
}

void ComputeDivevrgence(string input, TriMesh& mesh, UnitHole& hole)
{
#pragma region geodesic, controlVertex Settings
	mesh.update_normals();
	vector<vector<double>> geodesicDistance(mesh.n_faces());

	Calc_geodesicDistance(mesh, hole, geodesicDistance); // V, F에 데이터 입력

	// calc geodesicdistance

	VectorXd divMatri_x = VectorXd::Zero(mesh.n_vertices());
	VectorXd divMatri_y = VectorXd::Zero(mesh.n_vertices());
	VectorXd divMatri_z = VectorXd::Zero(mesh.n_vertices());

	ControlVertex CV;

	vector<bool> isBoundaryVertex;
	isBoundaryVertex.resize(mesh.n_vertices(), true);

	for (int i = 0; i < hole.new_Vertices.size(); i++)
	{
		isBoundaryVertex[hole.new_Vertices[i].idx()] = false;
	}

	CV.m_mesh = &mesh;

	for (int i = 0; i < hole.hole_original_vertices.size(); i++)
	{
		CV.m_HandVindex.push_back(hole.hole_original_vertices[i].first.idx());
		Point3D point(hole.hole_original_vertices[i].second[0],
			hole.hole_original_vertices[i].second[1],
			hole.hole_original_vertices[i].second[2]);
		CV.m_posOfHand.push_back(point);
	}

	//CV.ComputeFreeVertexWeight();

	PoissonDeformation PD(&mesh);	// poissonDeformation 객체 생성
	PD.poissonCtrl = CV;	// Original Vertices 를 PD 에서는 controlVertex로 대입
#pragma endregion

	int vid = 0, l = 0, r = 0;
	TriMesh::VertexHandle vh0, vh1, vh2;

	vector<pair<TriMesh::FaceHandle, Vector3D>> result_Fh_Vec3D; // face handle & result vec3d

	for (int i = 0; i < hole.new_Faces.size(); i++) // new Verteices
	{
		TriMesh::FaceVertexIter fv_it = mesh.fv_begin(hole.new_Faces[i]);
		vector<TriMesh::Point> face_vertices;
		vector<vec3> vertices_vec3;

		vid = mesh.fv_begin(hole.new_Faces[i])->idx();

		int tri_vid0 = fv_it.handle().idx(); fv_it++;
		int tri_vid1 = fv_it.handle().idx(); fv_it++;
		int tri_vid2 = fv_it.handle().idx();

		if (tri_vid0 == vid) { l = tri_vid1; r = tri_vid2; }
		else if (tri_vid1 == vid) { l = tri_vid0; r = tri_vid2; }
		else { l = tri_vid0; r = tri_vid1; }
		vh0 = mesh.vertex_handle(vid);
		vh1 = mesh.vertex_handle(l);
		vh2 = mesh.vertex_handle(r);

		face_vertices.push_back(mesh.point(vh0));
		face_vertices.push_back(mesh.point(vh1));
		face_vertices.push_back(mesh.point(vh2));

		vertices_vec3.push_back(vec3(face_vertices[0][0], face_vertices[0][1], face_vertices[0][2]));
		vertices_vec3.push_back(vec3(face_vertices[1][0], face_vertices[1][1], face_vertices[1][2]));
		vertices_vec3.push_back(vec3(face_vertices[2][0], face_vertices[2][1], face_vertices[2][2]));

		TriMesh::Normal n = mesh.normal(hole.new_Faces[i]);	// this interior triangle의 normal

		double x_n_average = 0;
		double y_n_average = 0;
		double z_n_average = 0;

		for (int j = 0; j < hole.original_Faces.size(); j++)
		{
			TriMesh::Normal FaceNormal = mesh.normal(hole.original_Faces[j]);

			x_n_average += (FaceNormal[0] / geodesicDistance[hole.new_Faces[i].idx()][j]);
			y_n_average += (FaceNormal[1] / geodesicDistance[hole.new_Faces[i].idx()][j]);
			z_n_average += (FaceNormal[2] / geodesicDistance[hole.new_Faces[i].idx()][j]);
		}

		x_n_average /= hole.original_Faces.size();
		y_n_average /= hole.original_Faces.size();
		z_n_average /= hole.original_Faces.size();

		TriMesh::Normal new_normal; // new normal
		new_normal[0] = x_n_average;
		new_normal[1] = y_n_average;
		new_normal[2] = z_n_average;

		vec3 f_n(n[0], n[1], n[2]);
		vec3 new_n(new_normal[0], new_normal[1], new_normal[2]);
		new_n.normalize();

		mat3 rotateAlign_mat3 = rotateAlign(f_n, new_n); // rotateAlign

		vec3 res = rotateAlign_mat3 * f_n;

		/*
		삼각형의 무게중심이 원점으로 가도록 이동하고 mat3행렬을 곱해서 회전하고 다시 원점이 원래 무게중심 위치로 가도록 이동
		*/

		vec3 G = get_G(face_vertices[0], face_vertices[1], face_vertices[2]);	// 무게 중심 계산

		vec3 result_vh0 = (rotateAlign_mat3 * (vertices_vec3[0] - G)) + G;
		vec3 result_vh1 = (rotateAlign_mat3 * (vertices_vec3[1] - G)) + G;
		vec3 result_vh2 = (rotateAlign_mat3 * (vertices_vec3[2] - G)) + G;


		//triangle local transform 
		Point3D source_, right_, left_;
		source_.m_x = result_vh0[0];
		source_.m_y = result_vh0[1];
		source_.m_z = result_vh0[2];

		left_.m_x = result_vh1[0];
		left_.m_y = result_vh1[1];
		left_.m_z = result_vh1[2];

		right_.m_x = result_vh2[0];
		right_.m_y = result_vh2[1];
		right_.m_z = result_vh2[2];

		//compute divergence
		Vector3D W = PD.ComputeTriangleDiv(source_, left_, right_);

		result_Fh_Vec3D.push_back(make_pair(hole.new_Faces[i], W));

		divMatri_x[vid] += W.m_x;
		divMatri_y[vid] += W.m_y;
		divMatri_z[vid] += W.m_z;

		W = PD.ComputeTriangleDiv(left_, right_, source_);

		divMatri_x[l] += W.m_x;
		divMatri_y[l] += W.m_y;
		divMatri_z[l] += W.m_z;

		W = PD.ComputeTriangleDiv(right_, source_, left_);

		divMatri_x[r] += W.m_x;
		divMatri_y[r] += W.m_y;
		divMatri_z[r] += W.m_z;
	}

	for (int i = 0; i < hole.hole_original_vertices.size(); i++)
	{
		TriMesh::Point p = mesh.point(hole.hole_original_vertices[i].first);

		divMatri_x[hole.hole_original_vertices[i].first.idx()] = p[0];
		divMatri_y[hole.hole_original_vertices[i].first.idx()] = p[1];
		divMatri_z[hole.hole_original_vertices[i].first.idx()] = p[2];
	}

	Deformation(mesh, hole, isBoundaryVertex, divMatri_x, divMatri_y, divMatri_z);
}


void Fill_Hole_Main(string inputPath, TriMesh& mesh)
{
	// boundary 찾기
	map<int, bool> isOriginal_Vertex;
	map<int, bool> isOriginal_Face;
	vector<pair<TriMesh::VertexHandle, TriMesh::Point>> holes;	// vertexhandle, points
	vector<pair<int, double>> index_angles; // point index, andgles;
	double length_average = 0;

	mesh.garbage_collection();

	statusSettings(mesh);

	eraseFace_hasTwoHalfEdge_boundary(mesh);

	statusSettings(mesh);

	set_isOriginal_Vertex_array(mesh, isOriginal_Vertex); // set Original Vertices
	set_isOriginal_Face_array(mesh, isOriginal_Face); // set Original Faces
	length_average = edge_length_average(mesh); // calc edges average length


	allquadrics::Quadric qfit;
	vector<UnitHole> Unit_holes = find_All_UnitHoles(mesh, false, qfit);

	int targethole_index;
	int Quadrics_index;
	int filled_hole_count = 0;
	save(mesh, "test01");

	while (true)
	{
		cout << "탐지된 구멍 개수 : " << Unit_holes.size() << endl;
		cout << "채운 구멍의 개수 : " << filled_hole_count << endl;
		cout << "메울 수 있는 구멍의 index :  ";

		for (int i = 0; i < Unit_holes.size(); i++)
		{
			if (!Unit_holes[i].is_filled)
				cout << Unit_holes[i].hole_number << " ";
		}

		cout << endl << "메울 구멍의 index를 입력하십시오 (메우지 않음 : -1) : ";

		cin >> targethole_index;

		if (targethole_index == -1)
			break;

		//==============================Quadrics=================================================

		vector<allquadrics::data_pnw> quadrics_hole_vertices_data;
		for (int i = 0; i < Unit_holes[targethole_index].hole_original_vertices.size(); i++)
		{
			vec3 point(Unit_holes[targethole_index].hole_original_vertices[i].second[0]
				, Unit_holes[targethole_index].hole_original_vertices[i].second[1]
				, Unit_holes[targethole_index].hole_original_vertices[i].second[2]);
			vec3 normal(Unit_holes[targethole_index].original_vertices_normals[i].second[0]
				, Unit_holes[targethole_index].original_vertices_normals[i].second[1]
				, Unit_holes[targethole_index].original_vertices_normals[i].second[2]);

			quadrics_hole_vertices_data.push_back({ point, normal, 1 });
		}

		vector<allquadrics::Quadric> qfits;
		allquadrics::Quadric q;

		//fitAllQuadricTypes(quadrics_hole_vertices_data, qfits);

		cout << endl << "quadrics 선택 (0 : general, 1 : rotationally symmetric, 2 : ) : ";

		cin >> Quadrics_index;

		bool useQuadrics = true;

		if (Quadrics_index == -1)
			useQuadrics = false;

		//==============================Quadrics=================================================

		fill_hole(mesh, Unit_holes[targethole_index].hole_original_vertices,
			Unit_holes[targethole_index].vertexIndex_angles, isOriginal_Vertex, length_average, q, useQuadrics);

		Unit_holes[targethole_index].new_Vertices = find_newVertices(mesh, Unit_holes[targethole_index], Unit_holes, isOriginal_Vertex);
		Unit_holes[targethole_index].new_Faces = find_newFaces(mesh, Unit_holes[targethole_index], Unit_holes, isOriginal_Face);

		find_unvalid_Vertices(mesh, Unit_holes[targethole_index]);


		Unit_holes[targethole_index].is_filled = true;
		filled_hole_count++;


		//=============================PoissonDeformation/0220===================================
		//Ax = B, b = PoissonDeformation::ComputeDivergence(), A = MeshLaplaceSolver::ComputeLalacianMatrixA()

		mesh.garbage_collection();
		save(mesh, "result_" + to_string(targethole_index) + "_fillhole", true); // save result

		ComputeDivevrgence(inputPath, mesh, Unit_holes[targethole_index]);

		save(mesh, "result_" + to_string(targethole_index) + "_MeshDeformation", true); // save result

		//=============================PoissonDeformation/0220===================================
		cout << "작업 완료" << endl << endl;
	}
}


//============================================================OpenMesh============================


void MyTriMesh::LoadMesh()
{

	CString m_filePath;
	CFileDialog file(true);
	file.m_ofn.lpstrFilter = _T("Mesh File(*.obj)\0*.obj\0Off File(*.off)\0*.off\0All File(*.*)\0*.*\0\0");
	file.m_ofn.lpstrTitle = _T("OPEN");
	if (file.DoModal() == IDOK)
	{
		m_filePath = file.GetPathName();
		OpenMesh::IO::read_mesh(m_mesh, m_filePath.GetBuffer(0));
	}

	init();
	GetAABBofObject();


	m_mesh.request_vertex_normals();
	m_mesh.request_face_normals();
	OpenMesh::IO::Options opt;

	// If the file did not provide vertex normals, then calculate them
	if (!opt.check(OpenMesh::IO::Options::VertexNormal) &&
		m_mesh.has_face_normals() && m_mesh.has_vertex_normals())
	{
		// let the mesh update the normals
		m_mesh.update_normals();
	}

	std::string filePath;
	filePath = m_filePath;

	igl::readOBJ(filePath, V, F);

	Eigen::VectorXi VS, FS, VT, FT;
	// The selected vertex is the source
	VS.resize(1);
	VS << 1;
	// All vertices are the targets
	VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
	Eigen::VectorXd d;
	//igl::exact_geodesic(V, F, VS, FS, VT, FT, d);
	//ReComputeVertexNormal();

	// Ax = B b = PoissonDeformation::ComputeDivergence(), A = MeshLaplaceSolver::ComputeLalacianMatrixA()
	//-------------------------0220 ----------------------------------------
	Fill_Hole_Main(filePath, m_mesh);
	//-------------------------0220 ----------------------------------------
}