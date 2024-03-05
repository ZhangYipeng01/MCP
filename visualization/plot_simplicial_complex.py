import plotly.graph_objects as go
import numpy as np

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

def process_data(file_path, threshold):
    # 初始化存储变量
    dslice1, dslice2, dslice3, dslice4 = [], [], [], []
    current_slice = None

    # 打开并读取文件
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # 检测当前行属于哪个slice
            if 'Slice' in line:
                # 移除冒号并获取slice编号
                slice_number = int(line.split(' ')[1].replace(':', ''))
                current_slice = slice_number
            else:
                # 解析数据和值
                data, value = eval(line)
                # 根据阈值和slice号存储数据
                if value <= threshold:
                    if current_slice == 1:
                        dslice1.append(data)
                    elif current_slice == 2:
                        dslice2.append(data)
                    elif current_slice == 3:
                        dslice3.append(data)
                    # 假设可能存在第四个slice，虽然示例中没有提到
                    elif current_slice == 4:
                        dslice4.append(data)
                # 其他条件根据需要添加

    # 返回处理后的数据
    return dslice1, dslice2, dslice3, dslice4

def visualize_simplicial_complex(vertices,edges,triangles,savename):
    # 提取顶点坐标
    x_coords = [v[1][0] for v in vertices]
    y_coords = [v[1][1] for v in vertices]
    z_coords = [0 for _ in vertices]  # Z坐标设为0

    # 三角形顶点索引
    i = [tri[0] for tri in triangles]
    j = [tri[1] for tri in triangles]
    k = [tri[2] for tri in triangles]

    # 创建 Plotly 图形
    fig = go.Figure()

    # 添加边
    for edge in edges:
        fig.add_trace(go.Scatter3d(
            x=[x_coords[edge[0]], x_coords[edge[1]]],
            y=[y_coords[edge[0]], y_coords[edge[1]]],
            z=[0, 0],
            mode='lines',
            line=dict(color='grey', width=5)))

    # 添加顶点
    fig.add_trace(go.Scatter3d(
        x=x_coords,
        y=y_coords,
        z=z_coords,
        mode='markers',
        marker=dict(size=4, color='black')))

    # 添加三角形
    fig.add_trace(go.Mesh3d(
        x=x_coords,
        y=y_coords,
        z=z_coords,
        i=i, j=j, k=k,
        color='lightpink', opacity=.25))

    # 设置场景，隐藏背景和坐标轴标签
    fig.update_layout(scene=dict(
        xaxis=dict(showbackground=False, showticklabels=False, title=''),
        yaxis=dict(showbackground=False, showticklabels=False, title=''),
        zaxis=dict(showbackground=False, showticklabels=False, title='')))

    # 显示图形
    fig.write_html(savename+".html")

def Combnumber(k,n): #Calculate the combination number

    numerator_index = 1
    denominator_index = n
    numerator = 1
    denominator = 1
    for i in range(k):
        numerator = numerator * numerator_index
        denominator = denominator * denominator_index
        numerator_index = numerator_index +1
        denominator_index = denominator_index -1
    return int(denominator/numerator)

def vertex_to_index(vertex):
    result = 0
    for i in range(len(vertex)):
        result += Combnumber(i+1,vertex[i])
    return result

#主程序

coords = TXTtoList('Terephthalic-acid.txt','float')
file_path = 'Terephthalic-acid_fslices.txt'  # 这需要替换为实际的文件路径
threshold = 2.25  # 根据需求设置阈值
dslice1, dslice2, dslice3, dslice4 = process_data(file_path, threshold)

vertices = []
edges = []
triangles = []
for simplex in dslice1:
    if len(simplex) == 1:
        for vertex in simplex:
            coordinate = [0.0,0.0]
            for point in vertex:
                coordinate = np.array(coordinate) + np.array(coords[point])
            vertices.append([vertex_to_index(vertex),coordinate.tolist()])
    if len(simplex) == 2:
        edge = []
        for vertex in simplex:
            edge.append(vertex_to_index(vertex))
        edges.append(edge)
    if len(simplex) == 3:
        triangle = []
        for vertex in simplex:
            triangle.append(vertex_to_index(vertex))
        triangles.append(triangle)

vertices_index_list = []
edges_get_index = []
triangles_get_index = []
for v in vertices:
    vertices_index_list.append(v[0])

for e in edges:
    index_v_0 = vertices_index_list.index(e[0])
    index_v_1 = vertices_index_list.index(e[1])
    edges_get_index.append([index_v_0,index_v_1])

for t in triangles:
    index_v_0 = vertices_index_list.index(t[0])
    index_v_1 = vertices_index_list.index(t[1])
    index_v_2 = vertices_index_list.index(t[2])
    triangles_get_index.append([index_v_0, index_v_1,index_v_2])

visualize_simplicial_complex(vertices,edges_get_index,triangles_get_index,'Terephthalic_Delaunay(k=1).txt')


vertices = []
edges = []
triangles = []
for simplex in dslice2:
    if len(simplex) == 1:
        for vertex in simplex:
            coordinate = [0.0,0.0]
            for point in vertex:
                coordinate = np.array(coordinate) + np.array(coords[point])
            vertices.append([vertex_to_index(vertex),coordinate.tolist()])
    if len(simplex) == 2:
        edge = []
        for vertex in simplex:
            edge.append(vertex_to_index(vertex))
        edges.append(edge)
    if len(simplex) == 3:
        triangle = []
        for vertex in simplex:
            triangle.append(vertex_to_index(vertex))
        triangles.append(triangle)

vertices_index_list = []
edges_get_index = []
triangles_get_index = []
for v in vertices:
    vertices_index_list.append(v[0])

for e in edges:
    index_v_0 = vertices_index_list.index(e[0])
    index_v_1 = vertices_index_list.index(e[1])
    edges_get_index.append([index_v_0,index_v_1])

for t in triangles:
    index_v_0 = vertices_index_list.index(t[0])
    index_v_1 = vertices_index_list.index(t[1])
    index_v_2 = vertices_index_list.index(t[2])
    triangles_get_index.append([index_v_0, index_v_1,index_v_2])

visualize_simplicial_complex(vertices,edges_get_index,triangles_get_index,'Terephthalic_Delaunay(k=2).txt')


vertices = []
edges = []
triangles = []
for simplex in dslice3:
    if len(simplex) == 1:
        for vertex in simplex:
            coordinate = [0.0,0.0]
            for point in vertex:
                coordinate = np.array(coordinate) + np.array(coords[point])
            vertices.append([vertex_to_index(vertex),coordinate.tolist()])
    if len(simplex) == 2:
        edge = []
        for vertex in simplex:
            edge.append(vertex_to_index(vertex))
        edges.append(edge)
    if len(simplex) == 3:
        triangle = []
        for vertex in simplex:
            triangle.append(vertex_to_index(vertex))
        triangles.append(triangle)

vertices_index_list = []
edges_get_index = []
triangles_get_index = []
for v in vertices:
    vertices_index_list.append(v[0])

for e in edges:
    index_v_0 = vertices_index_list.index(e[0])
    index_v_1 = vertices_index_list.index(e[1])
    edges_get_index.append([index_v_0,index_v_1])

for t in triangles:
    index_v_0 = vertices_index_list.index(t[0])
    index_v_1 = vertices_index_list.index(t[1])
    index_v_2 = vertices_index_list.index(t[2])
    triangles_get_index.append([index_v_0, index_v_1,index_v_2])

visualize_simplicial_complex(vertices,edges_get_index,triangles_get_index,'Terephthalic_Delaunay(k=3).txt')

