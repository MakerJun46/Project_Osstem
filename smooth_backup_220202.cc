#include <iostream>
#define no_init_all
#include <Windows.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
 // -------------------- OpenMesh

 //--------------------- allquadrics
#include "quadricfitting.h"
#include "view.h"
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


typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

#define M_PI 3.14159

struct UnitHole
{
	int hole_number;
	bool is_filled;
	vector<pair<MyMesh::VertexHandle, MyMesh::Normal>> original_vertices_normals;
	vector<pair<MyMesh::VertexHandle, MyMesh::Point>> hole_original_vertices;
	vector<MyMesh::VertexHandle> new_Vertices;
	vector<MyMesh::FaceHandle> new_Faces;
	vector<pair<int, double>> vertexIndex_angles;
};

void save(MyMesh mesh, string name, bool resultSave = false)
{
	try
	{
		if (resultSave)
		{
			if (!OpenMesh::IO::write_mesh(mesh, "Result/" + name + ".obj"))
			{
				std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			}
		}
		else
		{
			if (!OpenMesh::IO::write_mesh(mesh, "Log/" + name + ".obj"))
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
double edge_length_average(MyMesh& mesh)
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
void set_isOriginal_Vertex_array(MyMesh& mesh, map<int, bool>& isOriginal_Vertex)
{
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		isOriginal_Vertex[it->idx()] = true;
	}
}
void set_isOriginal_Face_array(MyMesh& mesh, map<int, bool>& isOriginal_Face)
{
	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		isOriginal_Face[it->idx()] = true;
	}
}

// find boundary vertices and sort loop
void boundary_vertex_Sort(MyMesh mesh, vector<MyMesh::VertexHandle> vhandle_boundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, vector<MyMesh::VertexHandle>& second_hole)
{
	holes.clear();

	if (vhandle_boundary.empty())
		return;

	holes.push_back(make_pair(vhandle_boundary[0], mesh.point(vhandle_boundary[0])));
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	MyMesh::HalfedgeHandle connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
	MyMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfEdge);
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
void boundary_Update(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, vector<MyMesh::VertexHandle>& second_hole)
{
	vector<MyMesh::VertexHandle> vhandle_boundary;

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

double calc_sector_angle_UseQaudrics(MyMesh& mesh, MyMesh::HalfedgeHandle _in_heh, vec3 normal)
{
	MyMesh::Normal v0, v1;
	MyMesh::Normal f_n(normal[0], normal[1], normal[2]);
	mesh.calc_sector_vectors(_in_heh, v0, v1);

	MyMesh::Scalar a = norm(v0);
	MyMesh::Scalar b = norm(v1);

	MyMesh::Scalar denom = norm(v0)*norm(v1);
	if (denom == MyMesh::Scalar(0))
	{
		return 0;
	}
	MyMesh::Scalar cos_a = dot(v0, v1) / denom;
	if (mesh.is_boundary(_in_heh))
	{
		MyMesh::Scalar sign_a = dot(cross(v0, v1), f_n);
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
double calc_angle(MyMesh& mesh, MyMesh::VertexHandle vh, bool useQuadrics, allquadrics::Quadric qfit)
{
	int idx = vh.idx();

	double deg;

	if (useQuadrics)
	{
		MyMesh::Point p = mesh.point(vh);

		vec3 point(p[0], p[1], p[2]);

		vec3 normal = partial_derivative(qfit, point);

		MyMesh::HalfedgeHandle prev_halfedge = mesh.prev_halfedge_handle(mesh.halfedge_handle(vh));

		deg = calc_sector_angle_UseQaudrics(mesh, prev_halfedge, normal);
	}
	else
	{
		MyMesh::HalfedgeHandle prev_halfedge = mesh.prev_halfedge_handle(mesh.halfedge_handle(vh));

		deg = mesh.calc_sector_angle(prev_halfedge) * 180 / M_PI;
	}

	if (deg < 0)
		deg += 360;

	return deg;
}
void angle_Update(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, vector<pair<int, double>>& index_angles, bool useQaudrics, allquadrics::Quadric qfit)
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
bool isClose_OtherVertex_Exist(vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, int a, int b)
{
	MyMesh::Point va1; // a 이전 vertex Point
	MyMesh::Point va2; // a 다음 vertex Point
	MyMesh::Point vb1; // b 이전 vertex Point
	MyMesh::Point vb2; // b 다음 vertex Point

	MyMesh::Point va = holes[a].second;
	MyMesh::Point vb = holes[b].second;

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
void fill_triangles_afterMerge(MyMesh& mesh, MyMesh::VertexHandle currentVh)
{
	MyMesh::HalfedgeHandle prev_halfedge = mesh.halfedge_handle(currentVh);

	vector<MyMesh::VertexHandle> TrianglePoints;

	int count = 0;

	while (count < 2)
	{
		prev_halfedge = mesh.prev_halfedge_handle(prev_halfedge);
		MyMesh::VertexHandle targetVh = mesh.from_vertex_handle(prev_halfedge);

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
void merge_close_vertex(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, map<int, bool> isOriginal_Vertex, double edge_length_average, int k, bool useQuadrics, allquadrics::Quadric qfit)
{
	for (int i = 0; i < holes.size(); i++)
	{
		vec3 p1(holes[i].second[0], holes[i].second[1], holes[i].second[2]);

		for (int j = 0; j < holes.size(); j++)
		{
			if ((i - j + holes.size()) % holes.size() <= 2 || (j - i + holes.size()) % holes.size() <= 2)
				continue;

			vec3 p2(holes[j].second[0], holes[j].second[1], holes[j].second[2]);

			double interval = sqrt(pow((p1[0] - p2[0]), 2) + pow((p1[1] - p2[1]), 2) + pow((p1[2] - p2[2]), 2));

			MyMesh::HalfedgeHandle hh1 = mesh.find_halfedge(holes[i].first, holes[j].first);
			MyMesh::HalfedgeHandle hh2 = mesh.find_halfedge(holes[j].first, holes[i].first);

			if (interval < edge_length_average && !isOriginal_Vertex[holes[j].first.idx()] && !isOriginal_Vertex[holes[i].first.idx()]
				&& hh1 == MyMesh::InvalidHalfedgeHandle && hh2 == MyMesh::InvalidHalfedgeHandle && !isClose_OtherVertex_Exist(holes, i, j))
			{
				vec3 center_coord((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2);

				if (useQuadrics)
					center_coord = calc_target_point(qfit, center_coord);

				MyMesh::Point p_new(center_coord[0], center_coord[1], center_coord[2]);


				mesh.set_point(holes[j].first, p_new);

				vector<vector<MyMesh::VertexHandle>> triangle_Points;

				for (auto it = mesh.vf_begin(holes[i].first); it != mesh.vf_end(holes[i].first); it++)	// merge vertices
				{
					MyMesh::FaceHandle fh = it;
					vector<MyMesh::VertexHandle> tmp;
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
				MyMesh::VertexHandle target = holes[i].first;
				MyMesh::VertexHandle merge_vh = holes[j].first;

				holes.erase(it);

				mesh.delete_vertex(target);

				for (auto vh : triangle_Points)
				{
					vector<MyMesh::VertexHandle> tmp;
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

				save(mesh, "after_merge_" + to_string(k));

				return;
			}
		}
	}
}

// fill hole
void fill_hole(MyMesh& mesh_findBoundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>> holes, vector<pair<int, double>> index_angles,
	map<int, bool> isOriginal_Vertex, double edge_length_average, allquadrics::Quadric qfit, bool useQuadrics)
{
	vector<vector<MyMesh::VertexHandle>> triangle_Points;
	vector<MyMesh::VertexHandle> second_hole;

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
		vector<MyMesh::VertexHandle> tmp;
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

			if (useQuadrics)
				new_coord = calc_target_point(qfit, new_coord);	// calc quadrics

			MyMesh::Point new_point(new_coord[0], new_coord[1], new_coord[2]);

			MyMesh::VertexHandle vh_tmp = mesh_findBoundary.add_vertex(new_point); // 새로운 점 생성

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

			if (useQuadrics)
			{
				a1 = calc_target_point(qfit, a1);	// calc quadrics
				b1 = calc_target_point(qfit, b1);
			}

			MyMesh::Point new_Point_a(a1[0], a1[1], a1[2]);
			MyMesh::Point new_Point_b(b1[0], b1[1], b1[2]);

			MyMesh::VertexHandle vh_tmp_a = mesh_findBoundary.add_vertex(new_Point_a);
			MyMesh::VertexHandle vh_tmp_b = mesh_findBoundary.add_vertex(new_Point_b);

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
			vector<MyMesh::VertexHandle> face_h;

			for (int j = 0; j < triangle_Points[i].size(); j++)
			{
				face_h.push_back(triangle_Points[i][j]);
			}

			mesh_findBoundary.add_face(face_h);
		}
		//--------------------------------------------------- vertex들을 이어서 면으로 만들기


		triangle_Points.clear();
		tmp.clear();

		save(mesh_findBoundary, to_string(k));
		if (newVertexCreated)
			merge_close_vertex(mesh_findBoundary, holes, isOriginal_Vertex, edge_length_average, k, useQuadrics, qfit);

		boundary_Update(mesh_findBoundary, holes, second_hole);
		angle_Update(mesh_findBoundary, holes, index_angles, useQuadrics, qfit);

		k++;
	}
}

// find new vertices, faces after fill hole process
bool is_other_hole_vertex(MyMesh::VertexHandle vh, int targethole_index, vector<UnitHole> unit_holes)
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
bool is_other_hole_face(MyMesh::FaceHandle fh, int targethole_index, vector<UnitHole> unit_holes)
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
vector<MyMesh::VertexHandle> find_newVertices(MyMesh& mesh, UnitHole hole, vector<UnitHole> unit_holes, map<int, bool> isOriginal_Vertex)
{
	vector<MyMesh::VertexHandle> newVertexHandles;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		if (isOriginal_Vertex[it->idx()] || is_other_hole_vertex(it, hole.hole_number, unit_holes))
			continue;
		else
			newVertexHandles.push_back(*it);
	}

	return newVertexHandles;
}
vector<MyMesh::FaceHandle> find_newFaces(MyMesh& mesh, UnitHole hole, vector<UnitHole> unit_holes, map<int, bool> isOriginal_Face)
{
	vector<MyMesh::FaceHandle> newFaceHandles;

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
vector<MyMesh::VertexHandle> Original_boundary_vertex_Sort(MyMesh& mesh, vector<MyMesh::VertexHandle>& vhandle_boundary)
{
	vector<MyMesh::VertexHandle> hole_vertices;

	hole_vertices.push_back(vhandle_boundary[0]);
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	MyMesh::HalfedgeHandle connected_halfedge = mesh.halfedge_handle(vertex_tmp);
	MyMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfedge);
	vhandle_boundary.erase(vhandle_boundary.begin());

	while (true)
	{
		auto target_it = find(vhandle_boundary.begin(), vhandle_boundary.end(), next_v);

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
vector<UnitHole> find_All_UnitHoles(MyMesh& mesh, bool useQuadrics, allquadrics::Quadric qfit)
{
	vector<UnitHole> vec_tmp;
	vector<MyMesh::VertexHandle> vhandle_boundary;
	vector<MyMesh::VertexHandle> vhandle_tmp;
	int hole_number = 0;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		if (mesh.is_boundary(it))
			vhandle_boundary.push_back(*it);
	}

	while (!vhandle_boundary.empty())
	{
		vhandle_tmp = Original_boundary_vertex_Sort(mesh, vhandle_boundary);
		UnitHole uh_tmp;

		for (auto vh : vhandle_tmp)
		{
			MyMesh::Normal vertex_normal = mesh.normal(vh);

			uh_tmp.original_vertices_normals.push_back(make_pair(vh, vertex_normal));
			uh_tmp.hole_original_vertices.push_back(make_pair(vh, mesh.point(vh)));
		}

		angle_Update(mesh, uh_tmp.hole_original_vertices, uh_tmp.vertexIndex_angles, useQuadrics, qfit);

		uh_tmp.hole_number = hole_number;
		uh_tmp.is_filled = false;

		vec_tmp.push_back(uh_tmp);

		hole_number++;
	}

	return vec_tmp;
}

void find_unvalid_Vertices(MyMesh& mesh, UnitHole& hole)
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

int main(int argc, char* argv[])
{
	// boundary 찾기

	MyMesh mesh_findBoundary;
	map<int, bool> isOriginal_Vertex;
	map<int, bool> isOriginal_Face;
	vector<pair<MyMesh::VertexHandle, MyMesh::Point>> holes;	// vertexhandle, points
	vector<pair<int, double>> index_angles; // point index, andgles;
	double length_average = 0;

	char* input = "../Smoothing/OstemTestCase/Crown_01_Hole.stl";

	if (!OpenMesh::IO::read_mesh(mesh_findBoundary, input))
	{
		std::cerr << "Error: Cannot read mesh from " << std::endl;
		return 1;
	}
	else
	{
		std::cout << "read mesh complete" << std::endl;
	}

	mesh_findBoundary.request_vertex_normals();
	set_isOriginal_Vertex_array(mesh_findBoundary, isOriginal_Vertex); // set Original Vertices
	set_isOriginal_Face_array(mesh_findBoundary, isOriginal_Face); // set Original Faces
	length_average = edge_length_average(mesh_findBoundary); // calc edges average length

	allquadrics::Quadric qfit;
	vector<UnitHole> Unit_holes = find_All_UnitHoles(mesh_findBoundary, false, qfit);

	int targethole_index;
	int Quadrics_index;
	int filled_hole_count = 0;

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
		allquadrics::fitAllQuadricTypes(quadrics_hole_vertices_data, qfits);

		cout << endl << "quadrics 선택 (0 : general, 1 : rotationally symmetric, 2 : ) : ";

		cin >> Quadrics_index;

		bool useQuadrics = true;

		if (Quadrics_index == -1)
			useQuadrics = false;

		//==============================Quadrics=================================================

		fill_hole(mesh_findBoundary, Unit_holes[targethole_index].hole_original_vertices,
			Unit_holes[targethole_index].vertexIndex_angles, isOriginal_Vertex, length_average, qfits[Quadrics_index], useQuadrics);

		Unit_holes[targethole_index].new_Vertices = find_newVertices(mesh_findBoundary, Unit_holes[targethole_index], Unit_holes, isOriginal_Vertex);
		Unit_holes[targethole_index].new_Faces = find_newFaces(mesh_findBoundary, Unit_holes[targethole_index], Unit_holes, isOriginal_Face);

		find_unvalid_Vertices(mesh_findBoundary, Unit_holes[targethole_index]);

		save(mesh_findBoundary, "result_" + to_string(targethole_index) + "_fillhole", true); // save result

		Unit_holes[targethole_index].is_filled = true;
		filled_hole_count++;

		cout << "작업 완료" << endl << endl;
	}

	return 0;
}