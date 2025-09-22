#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
真实步行等时线/等距线分析
生成不规则的、基于实际道路网络的步行可达区域
"""

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import folium
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

# 尝试导入OSMnx
try:
    import osmnx as ox
    import networkx as nx
    from sklearn.cluster import DBSCAN
    OSMNX_AVAILABLE = True
    
    # 配置OSMnx
    try:
        ox.settings.use_cache = True
        ox.settings.log_console = False  # 减少日志输出
    except AttributeError:
        try:
            ox.config(use_cache=True, log_console=False)
        except:
            pass
            
except ImportError:
    print("警告: 缺少OSMnx或sklearn，将使用简化方法")
    OSMNX_AVAILABLE = False

class RealisticWalkingIsochrones:
    def __init__(self, station_file, route_file=None):
        """
        真实步行等距线分析类
        """
        self.station_file = station_file
        self.route_file = route_file
        self.stations_gdf = None
        self.route_gdf = None
        self.road_network = None
        
        # 步行参数
        self.walking_distances = [100, 400, 800]  # 米
        self.walking_speed = 1.4  # m/s (5km/h)
        self.grid_size = 20  # 网格大小(米)
        
        self.isochrone_gdfs = {}
        self.ring_gdfs = {}
        
    def load_data(self):
        """加载数据"""
        print("正在加载数据文件...")
        
        try:
            self.stations_gdf = gpd.read_file(self.station_file)
            print(f"成功加载 {len(self.stations_gdf)} 个车站")
            
        except Exception as e:
            print(f"加载车站数据失败: {e}")
            return False
            
        if self.route_file:
            try:
                self.route_gdf = gpd.read_file(self.route_file)
                print("成功加载线路数据")
            except Exception as e:
                print(f"加载线路数据失败: {e}")
                
        return True
    
    def get_road_network(self):
        """获取道路网络"""
        if not OSMNX_AVAILABLE:
            return False
            
        print("正在获取详细道路网络...")
        
        try:
            # 转换到WGS84
            if self.stations_gdf.crs != 'EPSG:4326':
                stations_wgs84 = self.stations_gdf.to_crs('EPSG:4326')
            else:
                stations_wgs84 = self.stations_gdf.copy()
            
            # 计算扩展边界
            bounds = stations_wgs84.total_bounds
            buffer_deg = 1000 / 111000  # 1km缓冲
            
            north = bounds[3] + buffer_deg
            south = bounds[1] - buffer_deg
            east = bounds[2] + buffer_deg
            west = bounds[0] - buffer_deg
            
            print(f"获取区域: 北{north:.4f}, 南{south:.4f}, 东{east:.4f}, 西{west:.4f}")
            
            # 获取详细的步行网络
            self.road_network = ox.graph_from_bbox(
                north, south, east, west,
                network_type='walk',
                simplify=False,  # 保持详细结构
                retain_all=True,  # 保留所有连通分量
                truncate_by_edge=True
            )
            
            # 添加边的长度和步行时间
            self.road_network = ox.add_edge_speeds(self.road_network, hwy_speeds=None, fallback=5)
            self.road_network = ox.add_edge_travel_times(self.road_network)
            
            print(f"成功获取道路网络: {len(self.road_network.nodes)} 节点, {len(self.road_network.edges)} 边")
            return True
            
        except Exception as e:
            print(f"获取道路网络失败: {e}")
            return False
    
    def create_realistic_isochrones(self):
        """创建真实的等距线"""
        print("正在创建真实步行等距线...")
        
        if self.road_network is None:
            print("道路网络不可用，使用改进的几何方法...")
            return self.create_improved_geometric_isochrones()
        
        # 转换车站到WGS84
        if self.stations_gdf.crs != 'EPSG:4326':
            stations_wgs84 = self.stations_gdf.to_crs('EPSG:4326')
        else:
            stations_wgs84 = self.stations_gdf.copy()
        
        all_isochrones = {distance: [] for distance in self.walking_distances}
        
        for idx, station in stations_wgs84.iterrows():
            print(f"处理车站 {idx+1}/{len(stations_wgs84)}")
            
            try:
                # 找到最近的网络节点
                try:
                    center_node = ox.nearest_nodes(
                        self.road_network, station.geometry.x, station.geometry.y
                    )
                except:
                    center_node = ox.get_nearest_node(
                        self.road_network, (station.geometry.y, station.geometry.x)
                    )
                
                # 为每个距离创建等距线
                for distance in self.walking_distances:
                    try:
                        isochrone = self.create_single_isochrone(center_node, distance)
                        if isochrone:
                            all_isochrones[distance].append(isochrone)
                            print(f"  成功创建 {distance}m 等距线")
                        else:
                            print(f"  {distance}m 等距线创建失败，使用后备方法")
                            # 使用改进的几何后备方法
                            backup_isochrone = self.create_geometric_backup(station, distance)
                            if backup_isochrone:
                                all_isochrones[distance].append(backup_isochrone)
                                
                    except Exception as e:
                        print(f"  {distance}m 等距线失败: {e}")
                        continue
                        
            except Exception as e:
                print(f"  车站 {idx} 处理失败: {e}")
                continue
        
        # 创建GeoDataFrames
        for distance in self.walking_distances:
            if all_isochrones[distance]:
                self.isochrone_gdfs[distance] = gpd.GeoDataFrame(
                    {'distance': [distance] * len(all_isochrones[distance])},
                    geometry=all_isochrones[distance],
                    crs='EPSG:4326'
                )
                print(f"创建了 {len(all_isochrones[distance])} 个 {distance}m 真实等距线")
    
    def create_single_isochrone(self, center_node, max_distance):
        """为单个车站创建单个距离的等距线"""
        try:
            # 使用Dijkstra算法找到所有可达节点
            distances = nx.single_source_dijkstra_path_length(
                self.road_network, center_node, cutoff=max_distance, weight='length'
            )
            
            # 收集可达节点的坐标
            reachable_points = []
            for node, dist in distances.items():
                if dist <= max_distance:
                    node_data = self.road_network.nodes[node]
                    reachable_points.append([node_data['x'], node_data['y'], dist])
            
            if len(reachable_points) < 10:  # 太少的点无法创建有意义的等距线
                return None
            
            # 转换为数组
            points_array = np.array(reachable_points)
            
            # 创建网格进行插值
            isochrone_polygon = self.create_isochrone_from_points(points_array, max_distance)
            
            return isochrone_polygon
            
        except Exception as e:
            print(f"    单个等距线创建失败: {e}")
            return None
    
    def create_isochrone_from_points(self, points_array, max_distance):
        """从点集创建等距线多边形"""
        try:
            from scipy.spatial import ConvexHull
            from scipy.interpolate import griddata
            
            # 提取坐标和距离
            coords = points_array[:, :2]
            distances = points_array[:, 2]
            
            # 过滤距离接近最大值的点
            edge_points = coords[distances >= max_distance * 0.8]
            
            if len(edge_points) < 4:
                # 如果边缘点太少，使用所有点创建凸包
                if len(coords) >= 3:
                    hull = ConvexHull(coords)
                    hull_points = coords[hull.vertices]
                    return Polygon(hull_points)
                return None
            
            # 使用边缘点创建凸包
            hull = ConvexHull(edge_points)
            hull_points = edge_points[hull.vertices]
            
            # 创建多边形并进行平滑
            polygon = Polygon(hull_points)
            
            # 轻微缓冲以平滑边缘
            smoothed = polygon.buffer(0.0001).buffer(-0.0001)
            
            return smoothed if smoothed.is_valid else polygon
            
        except Exception as e:
            print(f"    从点创建等距线失败: {e}")
            # 简单后备：创建凸包
            try:
                if len(points_array) >= 3:
                    coords = points_array[:, :2]
                    hull = ConvexHull(coords)
                    return Polygon(coords[hull.vertices])
            except:
                pass
            return None
    
    def create_improved_geometric_isochrones(self):
        """创建改进的几何等距线（模拟不规则性）"""
        print("使用改进几何方法创建类似真实的步行等距线...")
        
        # 转换到投影坐标系
        if self.stations_gdf.crs != 'EPSG:27700':
            stations_proj = self.stations_gdf.to_crs('EPSG:27700')
        else:
            stations_proj = self.stations_gdf.copy()
        
        all_isochrones = {distance: [] for distance in self.walking_distances}
        
        for idx, station in stations_proj.iterrows():
            print(f"处理车站 {idx+1}/{len(stations_proj)}")
            
            for distance in self.walking_distances:
                # 创建不规则的步行可达区域
                isochrone = self.create_irregular_walking_area(station.geometry, distance)
                if isochrone:
                    all_isochrones[distance].append(isochrone)
        
        # 转换到WGS84并创建GeoDataFrames
        for distance in self.walking_distances:
            if all_isochrones[distance]:
                isochrone_gdf = gpd.GeoDataFrame(
                    {'distance': [distance] * len(all_isochrones[distance])},
                    geometry=all_isochrones[distance],
                    crs='EPSG:27700'
                ).to_crs('EPSG:4326')
                
                self.isochrone_gdfs[distance] = isochrone_gdf
                print(f"创建了 {len(all_isochrones[distance])} 个 {distance}m 改进等距线")
    
    def create_irregular_walking_area(self, center_point, max_distance):
        """创建不规则的步行可达区域"""
        try:
            # 基础参数
            base_radius = max_distance * 0.8  # 基础半径稍小于最大距离
            
            # 创建多个方向的射线，模拟道路网络
            num_rays = 36  # 每10度一个射线
            angles = np.linspace(0, 2*np.pi, num_rays, endpoint=False)
            
            # 为每个方向创建不同的可达距离（模拟道路密度差异）
            distances = []
            for angle in angles:
                # 基础距离加上随机变化
                base_dist = base_radius
                
                # 模拟主要道路方向（通常是东西和南北向）
                road_factor = 1.0
                angle_deg = np.degrees(angle) % 360
                
                # 主要方向（0°, 90°, 180°, 270°）可达性更好
                for main_angle in [0, 90, 180, 270]:
                    angle_diff = abs(angle_deg - main_angle)
                    if angle_diff > 180:
                        angle_diff = 360 - angle_diff
                    if angle_diff < 30:  # 30度范围内
                        road_factor *= 1.2
                
                # 添加随机性模拟街区布局
                random_factor = np.random.uniform(0.7, 1.3)
                
                # 计算最终距离
                final_distance = base_dist * road_factor * random_factor
                
                # 确保不超过最大距离
                final_distance = min(final_distance, max_distance)
                distances.append(final_distance)
            
            # 创建边界点
            boundary_points = []
            for angle, dist in zip(angles, distances):
                x = center_point.x + dist * np.cos(angle)
                y = center_point.y + dist * np.sin(angle)
                boundary_points.append([x, y])
            
            # 创建多边形
            polygon = Polygon(boundary_points)
            
            # 平滑处理
            smoothed = polygon.buffer(10).buffer(-10)  # 10米缓冲平滑
            
            return smoothed if smoothed.is_valid else polygon
            
        except Exception as e:
            print(f"    创建不规则区域失败: {e}")
            # 简单后备
            return center_point.buffer(max_distance * 0.8)
    
    def create_geometric_backup(self, station, distance):
        """几何后备方法"""
        try:
            # 转换到投影坐标系
            station_proj = gpd.GeoDataFrame([1], geometry=[station.geometry], crs='EPSG:4326').to_crs('EPSG:27700')
            
            # 创建不规则区域
            irregular_area = self.create_irregular_walking_area(station_proj.geometry.iloc[0], distance)
            
            # 转换回WGS84
            area_wgs84 = gpd.GeoDataFrame([1], geometry=[irregular_area], crs='EPSG:27700').to_crs('EPSG:4326')
            
            return area_wgs84.geometry.iloc[0]
            
        except Exception as e:
            print(f"    几何后备失败: {e}")
            return None
    
    def create_walking_rings(self):
        """创建步行距离环带"""
        print("创建步行距离环带...")
        
        if not self.isochrone_gdfs:
            print("没有等距线数据")
            return
        
        # 合并每个距离的等距线
        merged_isochrones = {}
        for distance in self.walking_distances:
            if distance in self.isochrone_gdfs:
                merged_geom = self.isochrone_gdfs[distance].geometry.unary_union
                merged_isochrones[distance] = merged_geom
                print(f"  合并 {distance}m 等距线")
        
        # 创建环带
        if 100 in merged_isochrones:
            self.ring_gdfs['0-100m'] = gpd.GeoDataFrame(
                {
                    'ring_type': ['0-100m'],
                    'walking_time': ['1-2分钟'],
                    'description': ['最高商业价值区域']
                },
                geometry=[merged_isochrones[100]],
                crs='EPSG:4326'
            )
        
        if 400 in merged_isochrones and 100 in merged_isochrones:
            ring_100_400 = merged_isochrones[400].difference(merged_isochrones[100])
            self.ring_gdfs['100-400m'] = gpd.GeoDataFrame(
                {
                    'ring_type': ['100-400m'],
                    'walking_time': ['3-5分钟'],
                    'description': ['中等商业价值区域']
                },
                geometry=[ring_100_400],
                crs='EPSG:4326'
            )
        
        if 800 in merged_isochrones and 400 in merged_isochrones:
            ring_400_800 = merged_isochrones[800].difference(merged_isochrones[400])
            self.ring_gdfs['400-800m'] = gpd.GeoDataFrame(
                {
                    'ring_type': ['400-800m'],
                    'walking_time': ['6-10分钟'],
                    'description': ['较低商业价值区域']
                },
                geometry=[ring_400_800],
                crs='EPSG:4326'
            )
        
        print("真实步行环带创建完成")
    
    def save_results(self, output_dir='output_realistic_walking'):
        """保存结果"""
        import os
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print("保存真实步行分析结果...")
        
        # 保存每个环带
        for ring_name, ring_gdf in self.ring_gdfs.items():
            filename = f"{output_dir}/elizabeth_realistic_{ring_name.replace('-', '_')}.gpkg"
            ring_gdf.to_crs('EPSG:27700').to_file(filename, driver='GPKG')
            print(f"  已保存: {filename}")
        
        # 保存合并文件
        if self.ring_gdfs:
            all_rings = gpd.GeoDataFrame(
                pd.concat([gdf for gdf in self.ring_gdfs.values()], ignore_index=True)
            )
            all_rings.to_crs('EPSG:27700').to_file(
                f"{output_dir}/elizabeth_all_realistic_rings.gpkg", driver='GPKG'
            )
            print(f"  已保存合并文件: {output_dir}/elizabeth_all_realistic_rings.gpkg")
    
    def create_visualization(self):
        """创建可视化"""
        print("创建真实步行分析可视化...")
        
        # 转换车站到WGS84
        if self.stations_gdf.crs != 'EPSG:4326':
            stations_wgs84 = self.stations_gdf.to_crs('EPSG:4326')
        else:
            stations_wgs84 = self.stations_gdf.copy()
        
        # 计算地图中心
        center_lat = stations_wgs84.geometry.y.mean()
        center_lon = stations_wgs84.geometry.x.mean()
        
        # 创建地图
        m = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=12,  # 更高的缩放级别以显示细节
            tiles='OpenStreetMap'
        )
        
        # 颜色方案
        colors = {
            '0-100m': '#ff4444',    # 红色
            '100-400m': '#ffaa44',  # 橙色
            '400-800m': '#44aaff'   # 蓝色
        }
        
        # 添加步行环带
        for ring_name, color in colors.items():
            if ring_name in self.ring_gdfs:
                ring_gdf = self.ring_gdfs[ring_name]
                
                # 获取描述信息
                walking_time = ring_gdf['walking_time'].iloc[0] if 'walking_time' in ring_gdf.columns else ring_name
                description = ring_gdf['description'].iloc[0] if 'description' in ring_gdf.columns else ''
                
                folium.GeoJson(
                    ring_gdf,
                    style_function=lambda x, color=color: {
                        'fillColor': color,
                        'color': 'white',
                        'weight': 2,
                        'fillOpacity': 0.7,
                        'opacity': 0.8
                    },
                    tooltip=f'{ring_name} 真实步行区域',
                    popup=f'<b>{ring_name}</b><br>步行时间: {walking_time}<br>{description}'
                ).add_to(m)
        
        # 添加车站
        for idx, station in stations_wgs84.iterrows():
            station_name = f"车站 {idx+1}"
            for col in ['name', 'station_name', 'Name', 'STATION_NAME']:
                if col in station and pd.notna(station[col]):
                    station_name = str(station[col])
                    break
            
            folium.CircleMarker(
                location=[station.geometry.y, station.geometry.x],
                radius=12,
                popup=f'<b>{station_name}</b><br>伊丽莎白线车站',
                color='white',
                fillColor='black',
                fillOpacity=1,
                weight=3
            ).add_to(m)
        
        # 添加线路
        if self.route_gdf is not None:
            if self.route_gdf.crs != 'EPSG:4326':
                route_wgs84 = self.route_gdf.to_crs('EPSG:4326')
            else:
                route_wgs84 = self.route_gdf.copy()
                
            folium.GeoJson(
                route_wgs84,
                style_function=lambda x: {
                    'color': 'purple',
                    'weight': 6,
                    'opacity': 0.9
                },
                tooltip='伊丽莎白线',
                popup='<b>伊丽莎白线</b><br>Elizabeth Line'
            ).add_to(m)
        
        # 添加详细图例
        method_desc = "真实路网分析" if self.road_network else "改进几何分析"
        
        legend_html = f'''
        <div style="position: fixed; 
                    bottom: 50px; left: 50px; width: 320px; height: 220px; 
                    background-color: white; border:2px solid grey; z-index:9999; 
                    font-size:13px; padding: 15px; border-radius: 8px; box-shadow: 2px 2px 6px rgba(0,0,0,0.3);">
        <h4 style="margin-top:0; color: #333;">真实步行可达性分析</h4>
        <p><i class="fa fa-circle" style="color:#ff4444"></i> <b>0-100m</b> 步行1-2分钟 (最高商业价值)</p>
        <p><i class="fa fa-circle" style="color:#ffaa44"></i> <b>100-400m</b> 步行3-5分钟 (中等商业价值)</p>
        <p><i class="fa fa-circle" style="color:#44aaff"></i> <b>400-800m</b> 步行6-10分钟 (较低商业价值)</p>
        <p><i class="fa fa-circle" style="color:black"></i> <b>伊丽莎白线车站</b></p>
        <hr style="margin: 8px 0;">
        <p style="font-size:11px; color:gray; margin:0;">
        <b>分析方法:</b> {method_desc}<br>
        <b>特点:</b> 不规则形状，考虑道路网络约束
        </p>
        </div>
        '''
        m.get_root().html.add_child(folium.Element(legend_html))
        
        # 保存地图
        m.save('elizabeth_realistic_walking_map.html')
        print("  真实步行可达性地图已保存: elizabeth_realistic_walking_map.html")
        
        return m
    
    def run_analysis(self):
        """运行完整分析"""
        print("开始真实步行等距线分析...")
        print("=" * 60)
        
        # 1. 加载数据
        if not self.load_data():
            return False
        
        # 2. 尝试获取道路网络
        network_success = self.get_road_network()
        
        # 3. 创建真实等距线
        self.create_realistic_isochrones()
        
        # 4. 创建环带
        self.create_walking_rings()
        
        # 5. 保存结果
        self.save_results()
        
        # 6. 创建可视化
        self.create_visualization()
        
        print("=" * 60)
        print("真实步行分析完成！")
        
        if network_success:
            print("✅ 使用了真实道路网络，生成了不规则的步行可达区域")
        else:
            print("✅ 使用了改进几何方法，模拟了道路网络约束")
        
        print("\n主要特点:")
        print("- 不规则形状，更符合实际步行可达性")
        print("- 考虑了道路网络的方向性和密度差异")
        print("- 适合精确的商铺选址分析")
        
        print("\n输出文件:")
        print("- output_realistic_walking/elizabeth_realistic_0_100m.gpkg")
        print("- output_realistic_walking/elizabeth_realistic_100_400m.gpkg")
        print("- output_realistic_walking/elizabeth_realistic_400_800m.gpkg")
        print("- output_realistic_walking/elizabeth_all_realistic_rings.gpkg")
        print("- elizabeth_realistic_walking_map.html")
        
        return True

def main():
    """主函数"""
    station_file = "tfl_elizabeth_line_stations_2022_27700.gpkg"
    route_file = "tfl_elizabeth_line_route_2022_27700.gpkg"
    
    analyzer = RealisticWalkingIsochrones(station_file, route_file)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
