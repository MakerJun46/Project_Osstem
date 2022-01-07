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

using namespace allquadrics;
using namespace std;

extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }

// ----------------------------------------------------------------------------

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

#define M_PI 3.14159

struct UnitHole
{
	int hole_number;
	vector<pair<MyMesh::VertexHandle, MyMesh::Normal>> original_vertices_normals;
	vector<pair<MyMesh::VertexHandle, MyMesh::Point>> hole_original_vertices;
	vector<pair<int, double>> vertexIndex_angles;
};

void save(MyMesh mesh, string name)
{
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "Log/" + name + ".obj"))
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
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

void boundary_vertex_Sort(MyMesh mesh, vector<MyMesh::VertexHandle> vhandle_boundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes)
{
	if (vhandle_boundary.empty())
		return;

	holes.clear();
	holes.push_back(make_pair(vhandle_boundary[0], mesh.point(vhandle_boundary[0])));
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	MyMesh::HalfedgeHandle connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
	MyMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfEdge);
	vhandle_boundary.erase(vhandle_boundary.begin());

	MyMesh::Point tmp3 = mesh.point(vertex_tmp);

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
					cout << "loop 를 확인하고 point를 정상적으로 정렬하였습니다." << endl;
				}
				else
				{
					std::cerr << "Error : loop가 정상적이지 않습니다. 시작점과 끝점이 연결되지 않습니다." << endl;
					break;
				}
				return;
			}
			else
			{
				auto p = mesh.point(vertex_tmp);

				cout << "Error : mesh에 포함되지 않은 Vertex가 있습니다. Vertex idx : " << vertex_tmp.idx() << endl;
				std::cerr << "Vertex Coordinate : " << p[0] << " " << p[1] << " " << p[2] << endl;
				break;
			}
		}

		auto target_it = find(vhandle_boundary.begin(), vhandle_boundary.end(), next_v);

		if (target_it == vhandle_boundary.end() && next_v == vertex_startPoint)
		{
			//cout << "구멍이 2개 이상임, 남은 vertex들은 다른 구멍의 vertex" << endl;
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
void boundary_Update(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes)
{
	vector<MyMesh::VertexHandle> vhandle_boundary;

	MyMesh::VertexIter v_it, v_end(mesh.vertices_end());

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		if (mesh.is_boundary(v_it))
		{
			vhandle_boundary.push_back(v_it);
		}
	}

	boundary_vertex_Sort(mesh, vhandle_boundary, holes);
}


vector<MyMesh::VertexHandle> Original_boundary_vertex_Sort(MyMesh& mesh, vector<MyMesh::VertexHandle>& vhandle_boundary)
{
	vector<MyMesh::VertexHandle> hole_vertices;

	hole_vertices.push_back(vhandle_boundary[0]);
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	MyMesh::HalfedgeHandle connected_halfedge = mesh.halfedge_handle(vertex_tmp);
	MyMesh::VertexHandle next_v = mesh.to_vertex_handle(connected_halfedge);
	vhandle_boundary.erase(vhandle_boundary.begin());

	MyMesh::Point tmp3 = mesh.point(vertex_tmp);

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


void boundary_Update_UnitHole(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes)
{
	vector<MyMesh::VertexHandle> vhandle_boundary;

	for (int i = 0; i < holes.size(); i++)
	{
		if (mesh.is_valid_handle(holes[i].first) && mesh.is_boundary(holes[i].first))
		{
			vhandle_boundary.push_back(holes[i].first);
		}
	}

	boundary_vertex_Sort(mesh, vhandle_boundary, holes);
}

double calc_angle(MyMesh& mesh, MyMesh::VertexHandle vh)
{
	int idx = vh.idx();

	MyMesh::HalfedgeHandle prev_halfedge = mesh.prev_halfedge_handle(mesh.halfedge_handle(vh));

	double deg = mesh.calc_sector_angle(prev_halfedge) * 180 / M_PI;

	if (deg < 0)
		deg += 360;

	return deg;
}
void angle_Update(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, vector<pair<int, double>>& index_angles)
{
	index_angles.clear();
	mesh.request_face_normals();

	for (int i = 0; i < holes.size(); i++)
	{
		index_angles.push_back(make_pair(i, calc_angle(mesh, holes[i].first)));
	}

	sort(index_angles.begin(), index_angles.end(), cmpSecond);
}


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
void merge_close_vertex(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, map<int, bool> isOriginal_Vertex, double edge_length_average, int k)
{
	for (int i = 0; i < holes.size(); i++)
	{
		vec3 p1(holes[i].second[0], holes[i].second[1], holes[i].second[2]);

		for (int j = 0; j < holes.size(); j++)
		{
			if (i == j || i - 1 == j || i + 1 == j || (i == holes.size() - 1 && j == 0))
				continue;

			vec3 p2(holes[j].second[0], holes[j].second[1], holes[j].second[2]);

			double interval = sqrt(pow((p1[0] - p2[0]), 2) + pow((p1[1] - p2[1]), 2) + pow((p1[2] - p2[2]), 2));

			MyMesh::HalfedgeHandle hh1 = mesh.find_halfedge(holes[i].first, holes[j].first);
			MyMesh::HalfedgeHandle hh2 = mesh.find_halfedge(holes[j].first, holes[i].first);

			if (interval < edge_length_average && !isOriginal_Vertex[holes[j].first.idx()] && !isOriginal_Vertex[holes[i].first.idx()]
				&& hh1 == MyMesh::InvalidHalfedgeHandle && hh2 == MyMesh::InvalidHalfedgeHandle && !isClose_OtherVertex_Exist(holes, i, j))
			{
				vec3 center_coord((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2);
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


void fill_holes(MyMesh& mesh_findBoundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>> holes, vector<pair<int, double>> index_angles,
	map<int, bool> isOriginal_Vertex, double edge_length_average, vector<UnitHole> Unit_Holes)
{
	vector<vector<MyMesh::VertexHandle>> triangle_Points;

	int k = 0;

	while (!holes.empty())
	{
		cout << endl;
		for (auto i : index_angles)
		{
			cout << i.first << "번 사이각 : " << i.second << endl;
		}

		int index = index_angles[0].first;	// 사이각이 가장 작은 개체(내림차순 정렬의 0번째)
		vector<MyMesh::VertexHandle> tmp;
		vec3 index_vec(holes[index].second[0], holes[index].second[1], holes[index].second[2]);

		if (holes.size() < 3)
			break;

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
		merge_close_vertex(mesh_findBoundary, holes, isOriginal_Vertex, edge_length_average, k);
		//boundary_Update(mesh_findBoundary, holes);
		boundary_Update_UnitHole(mesh_findBoundary, holes);
		angle_Update(mesh_findBoundary, holes, index_angles);

		k++;
	}
}


vector<MyMesh::VertexHandle> find_newVertices(MyMesh& mesh, map<int, bool> isOriginal_Vertex)
{
	vector<MyMesh::VertexHandle> newVertexHandles;

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		if (isOriginal_Vertex[it->idx()])
			continue;
		else
			newVertexHandles.push_back(*it);
	}

	return newVertexHandles;
}
vector<MyMesh::FaceHandle> find_newFaces(MyMesh& mesh, map<int, bool> isOriginal_Face)
{
	vector<MyMesh::FaceHandle> newFaceHandles;

	for (auto it = mesh.faces_begin(); it != mesh.faces_end(); it++)
	{
		if (isOriginal_Face[it->idx()])
			continue;
		else
			newFaceHandles.push_back(*it);
	}

	return newFaceHandles;
}



vector<UnitHole> find_All_UnitHoles(MyMesh& mesh)
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
		cout << vhandle_boundary.size() << endl;
		vhandle_tmp = Original_boundary_vertex_Sort(mesh, vhandle_boundary);
		UnitHole uh_tmp;

		for (auto vh : vhandle_tmp)
		{
			MyMesh::Normal vertex_normal = mesh.normal(vh);

			uh_tmp.original_vertices_normals.push_back(make_pair(vh, vertex_normal));
			uh_tmp.hole_original_vertices.push_back(make_pair(vh, mesh.point(vh)));
		}

		angle_Update(mesh, uh_tmp.hole_original_vertices, uh_tmp.vertexIndex_angles);

		uh_tmp.hole_number = hole_number;

		vec_tmp.push_back(uh_tmp);

		hole_number++;
	}

	return vec_tmp;
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


	if (!OpenMesh::IO::read_mesh(mesh_findBoundary, "../Smoothing/bunny_multiHoles.obj"))
	{
		std::cerr << "Error: Cannot read mesh from " << std::endl;
		return 1;
	}
	else
	{
		std::cout << "read mesh complete" << std::endl;
	}

	mesh_findBoundary.request_vertex_normals();
	vector<UnitHole> Unit_holes = find_All_UnitHoles(mesh_findBoundary);

	set_isOriginal_Vertex_array(mesh_findBoundary, isOriginal_Vertex); // set Original Vertices
	//set_isOriginal_Face_array(mesh_findBoundary, isOriginal_Face); // set Original Faces
	length_average = edge_length_average(mesh_findBoundary); // calc edges average length

	//boundary_Update(mesh_findBoundary, holes);  // find boundary vertices and sort loop

	//angle_Update(mesh_findBoundary, holes, index_angles); // calc angles and sort
	//fill_holes(mesh_findBoundary, holes, index_angles, isOriginal_Vertex, length_average); // fill holes

	//vector<MyMesh::VertexHandle> newVertexHandles = find_newVertices(mesh_findBoundary, isOriginal_Vertex); // set new Vertices after fill holes
	//vector<MyMesh::FaceHandle> newFaceHandles = find_newFaces(mesh_findBoundary, isOriginal_Face); // set new Faces after fill holes

	fill_holes(mesh_findBoundary, Unit_holes[0].hole_original_vertices, Unit_holes[0].vertexIndex_angles, isOriginal_Vertex, length_average, Unit_holes);


	save(mesh_findBoundary, "result"); // save result

	cout << "작업 완료" << endl;

	return 0;
}