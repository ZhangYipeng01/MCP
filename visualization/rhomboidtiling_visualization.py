import re
import numpy as np
from itertools import combinations
import plotly.graph_objects as go
from scipy.spatial import ConvexHull
import ast

def TXTtoList(filepath, return_type='str'):
    fp = open(filepath, 'r', encoding='utf-8')
    content = fp.read()
    fp.close()
    rowlist = content.splitlines()
    recordlist = [row.split() for row in rowlist if row != 'END' if row.strip()]
    if return_type == 'str':
        return recordlist
    if return_type == 'int':
        int_recordlist = [[int(element) for element in inner_list] for inner_list in recordlist]
        return int_recordlist
    if return_type == 'float':
        float_recordlist = [[float(element) for element in inner_list] for inner_list in recordlist]
        return float_recordlist

def vertex_get_coords(vertex):
    if len(points_set[0]) == 2:
        coords = np.array([0.0, 0.0, 0.0])
    if len(points_set[0]) == 3:
        coords = np.array([0.0, 0.0, 0.0, 0.0])
    for index in vertex:
        coords = np.array(coords) + np.array(points_set[index]+[-1.0])
    return coords.tolist()

def plot_3d_convex_hull(points, labels, edges):
    hull = ConvexHull(points)
    fig = go.Figure()

    # 添加 Convex Hull 的表面
    fig.add_trace(go.Mesh3d(
        x=points[:,0], y=points[:,1], z=points[:,2],
        i=hull.simplices[:,0], j=hull.simplices[:,1], k=hull.simplices[:,2],
        opacity=0.5, color='lightblue'
    ))

    # 绘制所有边
    for edge in edges:
        fig.add_trace(go.Scatter3d(
            x=[edge[0][0], edge[1][0]],
            y=[edge[0][1], edge[1][1]],
            z=[edge[0][2], edge[1][2]],
            mode='lines', line=dict(color='red', width=2)
        ))

    # 添加顶点和顶点标签
    fig.add_trace(go.Scatter3d(
        x=points[:, 0], y=points[:, 1], z=points[:, 2],
        mode='markers+text',
        marker=dict(size=5, color='black'),
        text=labels,
        textposition="top center",
        textfont=dict(
            size=20,
            color="black"
        )
    ))

    # 隐藏坐标轴和背景，隐藏坐标轴标签
    fig.update_layout(
        showlegend=False,  # 隐藏图例
        scene=dict(
            xaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏X轴标签
            yaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏Y轴标签
            zaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏Z轴标签
        ),
        width=700,
        margin=dict(r=20, l=10, b=10, t=10)
    )

    fig.show()

def plot_rhomboidtiling(rhomboids, all_points, all_labels, edges):
    fig = go.Figure()

    # 添加 Rhomboid tiling 的表面
    for points in rhomboids:
        hull = ConvexHull(points)
        nppoints = np.array(points)
        fig.add_trace(go.Mesh3d(
            x=nppoints[:,0], y=nppoints[:,1], z=nppoints[:,2],
            i=hull.simplices[:,0], j=hull.simplices[:,1], k=hull.simplices[:,2],
            opacity=0.5, color='lightblue'
        ))

    # 绘制所有边
    for edge in edges:
        fig.add_trace(go.Scatter3d(
            x=[edge[0][0], edge[1][0]],
            y=[edge[0][1], edge[1][1]],
            z=[edge[0][2], edge[1][2]],
            mode='lines', line=dict(color='red', width=2)
        ))

    # 添加顶点和顶点标签
    fig.add_trace(go.Scatter3d(
        x=all_points[:, 0], y=all_points[:, 1], z=all_points[:, 2],
        mode='markers+text',
        marker=dict(size=5, color='black'),
        text=all_labels,
        textposition="top center",
        textfont=dict(
            size=20,
            color="black"
        )
    ))

    # 隐藏坐标轴和背景，隐藏坐标轴标签
    fig.update_layout(
        showlegend=False,  # 隐藏图例
        scene=dict(
            xaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏X轴标签
            yaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏Y轴标签
            zaxis=dict(showbackground=False, showticklabels=False, zeroline=False, title=''),  # 隐藏Z轴标签
        ),
        width=700,
        margin=dict(r=20, l=10, b=10, t=10)
    )

    fig.show()



points_set = TXTtoList('4points.txt','float')
print(points_set)

# 初始化一个空列表来存储edges_core
edges_core = []


with open('4points_rhomboidtiling.txt', 'r') as file:
    for line in file:
        if 'd:1' in line:
            # 使用正则表达式匹配 vxs: 后的列表数据
            match = re.search(r'vxs:\s*(\[\[.*?\]\])', line)
            if match:
                # 获取匹配到的列表字符串
                vxs_data = match.group(1)
                # 由于我们不再使用 eval(), 这里可以采用一个简单的方法来解析这个字符串
                # 注意：这种方法假设数据格式非常规范，仅适用于格式良好的数据
                # 对于更复杂或不规范的数据格式，可能需要更复杂的解析方法
                try:
                    edges = ast.literal_eval(vxs_data)
                    edges_core.append(edges)
                except ValueError as e:
                    print(f"Error parsing vxs_data: {vxs_data}")
                    print(e)

edges_coord_set = []
for point in points_set:
    edges_coord_set.append([[0,0,0],point+[-1]])
for edge_indexes in edges_core:
    v_0 = vertex_get_coords(edge_indexes[0])
    v_1 = vertex_get_coords(edge_indexes[1])
    edges_coord_set.append([v_0,v_1])

rhomboid_core = []
# 打开并读取文件
with open('4points_rhomboids.txt', 'r') as file:
    for line in file:
        line = line.strip()
        lists = re.findall(r'\[(.*?)\]', line)

        temp_result = []
        for list_str in lists:
            if list_str:
                list_int = [int(x) for x in list_str.split(',')]
            else:
                list_int = []
            temp_result.append(list_int)
        rhomboid_core.append(temp_result)

rhomboids_tiling_vertexes_set = []
all_vertexes = []
labels_set = []
all_labels_set = []
for rhomboid in rhomboid_core:
    labels = []
    rhomboid_vertexes_set = []
    X_in = rhomboid[0]
    X_on = rhomboid[1]
    v_X_in = vertex_get_coords(X_in)
    if len(X_in) == 0:
        labels.append('O')
    else:
        labels.append(str(X_in))
    rhomboid_vertexes_set.append(v_X_in)
    if v_X_in not in all_vertexes:
        if len(X_in) == 0:
            all_labels_set.append('O')
        else:
            all_labels_set.append(str(X_in))
        all_vertexes.append(v_X_in)

    for r in range(len(X_on)):
        for subset in combinations(X_on, r+1):
            if len(X_in) == 0:
                labels.append(str(list(subset)))
            else:
                labels.append(str(X_in+list(subset)))
            root_point = np.array(v_X_in) + np.array(vertex_get_coords(subset))
            rhomboid_vertexes_set.append(root_point.tolist())
            if root_point.tolist() not in all_vertexes:
                if len(X_in) == 0:
                    all_labels_set.append(str(list(subset)))
                else:
                    all_labels_set.append(str(X_in+list(subset)))
                all_vertexes.append(root_point.tolist())

    rhomboids_tiling_vertexes_set.append(rhomboid_vertexes_set)
    labels_set.append(labels)

for i in range(len(rhomboids_tiling_vertexes_set)):
    plot_3d_convex_hull(np.array(rhomboids_tiling_vertexes_set[i]),labels_set[i],edges_coord_set)

plot_rhomboidtiling(rhomboids_tiling_vertexes_set,np.array(all_vertexes),all_labels_set,edges_coord_set)
#fig.show()

