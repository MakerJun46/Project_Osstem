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
	vector<pair<MyMesh::VertexHandle, MyMesh::Point>> hole_vertices;
	vector<pair<int, double>> index_angles;
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

void set_isOriginal_Vertex_array(MyMesh& mesh, map<int, bool>& isOriginal_Vertex)
{
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
	{
		isOriginal_Vertex[it->idx()] = true;
	}
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

double angle(MyMesh& mesh, MyMesh::VertexHandle vh, MyMesh::Point p1, MyMesh::Point p_target, MyMesh::Point p2)
{
	vec3 a(p1[0], p1[1], p1[2]);
	vec3 x(p_target[0], p_target[1], p_target[2]);
	vec3 b(p2[0], p2[1], p2[2]);

	a -= x;
	b -= x;

	a.normalize();
	b.normalize();

	vec3 cross_Product(a[1] * b[2] - a[2] * a[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0]);

	auto it = mesh.vf_begin(vh);

	MyMesh::Normal face_normal = mesh.normal(it);

	double inner_Product = face_normal[0] * cross_Product[0] +
		face_normal[1] * cross_Product[1] +
		face_normal[2] * cross_Product[2];


	double rad = acos(a * b);

	double deg = rad * 180 / M_PI;

	if (inner_Product < 0)
		deg = 360 - deg;

	return deg;
}

void angle_Update(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, vector<pair<int, double>>& index_angles)
{
	index_angles.clear();
	mesh.request_face_normals();
	for (int i = 0; i < holes.size(); i++)
	{
		if (i == 0)
		{
			index_angles.push_back(make_pair(i, angle(mesh, holes[i].first, holes[holes.size() - 1].second, holes[i].second, holes[i + 1].second)));
		}
		else if (i == holes.size() - 1)
		{
			index_angles.push_back(make_pair(i, angle(mesh, holes[i].first, holes[i - 1].second, holes[i].second, holes[0].second)));
		}
		else
		{
			index_angles.push_back(make_pair(i, angle(mesh, holes[i].first, holes[i - 1].second, holes[i].second, holes[i + 1].second)));
		}
	}

	sort(index_angles.begin(), index_angles.end(), cmpSecond);
}

void vertex_Sort(MyMesh mesh, vector<MyMesh::VertexHandle> vhandle_boundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes)
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

	vector<UnitHole> unit_holes;
	UnitHole uh_tmp;

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

				uh_tmp.hole_vertices = holes;
				unit_holes.push_back(uh_tmp);

				if (next_v.idx() == vertex_startPoint.idx())
				{
					cout << "loop ?? ???????? point?? ?????????? ??????????????." << endl;
					cout << "???? holes ???? : " << unit_holes.size() << "??" << endl;
				}
				else
				{
					std::cerr << "Error : loop?? ?????????? ????????. ???????? ?????? ???????? ????????." << endl;
					break;
				}
				return;
			}
			else
			{
				auto p = mesh.point(vertex_tmp);

				cout << "Error : mesh?? ???????? ???? Vertex?? ????????. Vertex idx : " << vertex_tmp.idx() << endl;
				std::cerr << "Vertex Coordinate : " << p[0] << " " << p[1] << " " << p[2] << endl;
				break;
			}
		}

		auto target_it = find(vhandle_boundary.begin(), vhandle_boundary.end(), next_v);

		if (target_it == vhandle_boundary.end() && next_v == vertex_startPoint)
		{
			cout << "?????? 2?? ??????, ???? vertex???? ???? ?????? vertex" << endl;

			uh_tmp.hole_vertices = holes;
			unit_holes.push_back(uh_tmp);

			holes.clear();
			holes.push_back(make_pair(vhandle_boundary[0], mesh.point(vhandle_boundary[0])));
			vertex_startPoint = vhandle_boundary[0];
			vertex_tmp = vhandle_boundary[0];
			vhandle_boundary.erase(vhandle_boundary.begin());
			connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
			next_v = mesh.to_vertex_handle(connected_halfEdge);
		}
		else
		{
			holes.push_back(make_pair(*target_it, mesh.point(*target_it)));
			vertex_tmp = *target_it;

			connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
			next_v = mesh.to_vertex_handle(connected_halfEdge);

			vhandle_boundary.erase(target_it);
		}


		//for (auto it = vhandle_boundary.begin(); it != vhandle_boundary.end(); it++)
		//{
		//	if (next_v.idx() == it->idx())
		//	{
		//		holes.push_back(make_pair(*it, mesh.point(*it)));
		//		vertex_tmp = *it;

		//		connected_halfEdge = mesh.halfedge_handle(vertex_tmp);
		//		next_v = mesh.to_vertex_handle(connected_halfEdge);

		//		vhandle_boundary.erase(it);
		//		break;
		//	}
		//}
	}
}

bool isClose_OtherVertex_Exist(vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, int a, int b)
{
	MyMesh::Point va1; // a ???? vertex Point
	MyMesh::Point va2; // a ???? vertex Point
	MyMesh::Point vb1; // b ???? vertex Point
	MyMesh::Point vb2; // b ???? vertex Point

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

	vertex_Sort(mesh, vhandle_boundary, holes);
}

void merge_close_vertex(MyMesh& mesh, vector<pair<MyMesh::VertexHandle, MyMesh::Point>>& holes, map<int, bool> isOriginal_Vertex, double edge_length_average)
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

				mesh.delete_vertex(holes[i].first);

				mesh.garbage_collection();

				for (auto vh : triangle_Points)
				{
					vector<MyMesh::VertexHandle> tmp;
					for (auto v : vh)
					{
						tmp.push_back(v);
					}

					cout << "p1" << endl;
					mesh.add_face(tmp);
					tmp.clear();
				}

				return;
			}
		}
	}
}

void fill_holes(MyMesh& mesh_findBoundary, vector<pair<MyMesh::VertexHandle, MyMesh::Point>> holes, vector<pair<int, double>> index_angles,
	map<int, bool> isOriginal_Vertex, double edge_length_average)
{
	vector<vector<MyMesh::VertexHandle>> triangle_Points;

	int k = 0;

	while (!holes.empty())
	{
		cout << endl;
		for (auto i : index_angles)
		{
			cout << i.first << "?? ?????? : " << i.second << endl;
		}

		int index = index_angles[0].first;	// ???????? ???? ???? ????(???????? ?????? 0????)
		vector<MyMesh::VertexHandle> tmp;
		vec3 index_vec(holes[index].second[0], holes[index].second[1], holes[index].second[2]);

		if (holes.size() < 3)
			break;

		if (index_angles[0].second <= 75.1) // 75???? ???? ???? ???? -> ?? 0?? ????, ?? 1?? ????, ?????? 1?? ????
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

			triangle_Points.push_back(tmp); // ?????? ????

			holes.erase(find(holes.begin(), holes.end(), holes[index]));
		}

		else if (75.1 < index_angles[0].second && index_angles[0].second <= 135.1) // 75 < ?????? < 135 ?? ???? -> ?? 1?? ????, ?? 1?? ????, ?????? 2?? ????
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
			vec3 new_vec = v1 + v2; // ?????? ?? ???? (a ???? + b ????)

			new_vec.normalize();

			new_vec = new_vec;// *((v1.length() + v2.length()) / 2);

			new_vec *= max(edge_length_average, ((v1.length() + v2.length()) / 2));

			vec3 new_coord = index_vec + new_vec;

			MyMesh::Point new_point(new_coord[0], new_coord[1], new_coord[2]);

			MyMesh::VertexHandle vh_tmp = mesh_findBoundary.add_vertex(new_point); // ?????? ?? ????

			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp);
			triangle_Points.push_back(tmp); // ?????? 1 ????

			tmp.clear();
			tmp.push_back(vh_tmp);
			tmp.push_back(holes[index].first);

			if (index == holes.size() - 1)
				tmp.push_back(holes[0].first);
			else
				tmp.push_back(holes[index + 1].first);

			triangle_Points.push_back(tmp); // ?????? 2 ????

			holes[index] = make_pair(vh_tmp, new_point);// ?????? ?????????? holes?? angle?? ???? ?? ?????? ??????
		}

		else // 135 < ?????? ?? ???? -> ?? 2?? ????, ?? 1?? ????, ?????? 3?? ????
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

			center_coord = center_coord * ((v1.length() + v2.length()) / 2); //center coord = v1?? v2?? ????


			vec3 a1 = (v1 + center_coord * 2).normalize();// *(v1.length() + center_coord.length() * 2) / 3;
			vec3 b1 = (center_coord * 2 + v2).normalize();// *(center_coord.length() * 2 + v2.length()) / 3;
			a1 *= max(edge_length_average, (v1.length() + center_coord.length() * 2) / 3);
			b1 *= max(edge_length_average, (center_coord.length() * 2 + v2.length()) / 3);

			/*
			vec3 a1 = v1 + center_coord;	// v1?? center coord?? ????
			vec3 b1 = v2 + center_coord;	// v2?? center coord?? ????

			a1.normalize();
			b1.normalize();

			a1 = a1 * ((v1.length() + center_coord.length()) / 2);
			b1 = b1 * ((v2.length() + center_coord.length()) / 2);
			*/


			a1 += index_vec;
			b1 += index_vec;

			MyMesh::Point new_Point_a(a1[0], a1[1], a1[2]);
			MyMesh::Point new_Point_b(b1[0], b1[1], b1[2]);

			MyMesh::VertexHandle vh_tmp_a = mesh_findBoundary.add_vertex(new_Point_a);
			MyMesh::VertexHandle vh_tmp_b = mesh_findBoundary.add_vertex(new_Point_b);

			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp_a);
			triangle_Points.push_back(tmp);	// ?????? 1 ????

			tmp.clear();
			tmp.push_back(vh_tmp_a);
			tmp.push_back(holes[index].first);
			tmp.push_back(vh_tmp_b);
			triangle_Points.push_back(tmp); // ?????? 2 ????

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

			triangle_Points.push_back(tmp);	// ?????? 3 ????

			holes.insert(holes.begin() + index, make_pair(vh_tmp_b, new_Point_b)); // ?????? point2 ????
			holes.insert(holes.begin() + index, make_pair(vh_tmp_a, new_Point_a));	// ?????? point1 ????
			holes.erase(holes.begin() + index + 2); // point ????
		}

		//--------------------------------------------------- vertex???? ?????? ?????? ??????
		for (int i = 0; i < triangle_Points.size(); i++)
		{
			vector<MyMesh::VertexHandle> face_h;

			for (int j = 0; j < triangle_Points[i].size(); j++)
			{
				face_h.push_back(triangle_Points[i][j]);
			}

			cout << "p2" << endl;
			mesh_findBoundary.add_face(face_h);
		}
		//--------------------------------------------------- vertex???? ?????? ?????? ??????

		triangle_Points.clear();
		tmp.clear();

		merge_close_vertex(mesh_findBoundary, holes, isOriginal_Vertex, edge_length_average);
		save(mesh_findBoundary, to_string(k));
		boundary_Update(mesh_findBoundary, holes);
		angle_Update(mesh_findBoundary, holes, index_angles);

		k++;
	}
}



int main(int argc, char* argv[])
{
	// boundary ????

	MyMesh mesh_findBoundary;
	map<int, bool> isOriginal_Vertex;	// ?????? vertex???? ????
	vector<pair<MyMesh::VertexHandle, MyMesh::Point>> holes;	// vertexhandle, points
	vector<pair<int, double>> index_angles; // point index, andgles;
	double length_average = 0;


	if (!OpenMesh::IO::read_mesh(mesh_findBoundary, "../Smoothing/testObj.obj"))
	{
		std::cerr << "Error: Cannot read mesh from " << std::endl;
		return 1;
	}
	else
	{
		std::cout << "read mesh complete" << std::endl;
	}

	mesh_findBoundary.request_face_colors();
	mesh_findBoundary.request_vertex_colors();
	mesh_findBoundary.request_face_normals(); // ?????? ???? request


	set_isOriginal_Vertex_array(mesh_findBoundary, isOriginal_Vertex); // ?????? vertex???? index ????
	length_average = edge_length_average(mesh_findBoundary); // edge?? ???? ???? ????

	boundary_Update(mesh_findBoundary, holes);  // boundary ???? sort
	angle_Update(mesh_findBoundary, holes, index_angles); // ?????? ???? ?? sort
	fill_holes(mesh_findBoundary, holes, index_angles, isOriginal_Vertex, length_average); // ???? ?????? ????

	save(mesh_findBoundary, "result"); // ???? ????

	cout << "???? ????" << endl;

	return 0;
}