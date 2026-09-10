import sys
import os
import re
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QSlider, 
                               QFileDialog, QGroupBox, QFormLayout, QLineEdit, 
                               QMessageBox, QGraphicsRectItem, QComboBox, QScrollArea,
                               QDialog, QTableWidget, QTableWidgetItem, QHeaderView,
                               QPlainTextEdit, QSplitter)
from PySide6.QtCore import Qt, QProcess, QTimer
import pyqtgraph as pg
import pyvista as pv
import numpy as np

class PreferencesDialog(QDialog):
    def __init__(self, parent, current_size, current_theme):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setFixedWidth(350)
        layout = QVBoxLayout(self)
        
        form = QFormLayout()
        self.font_slider = QSlider(Qt.Horizontal)
        self.font_slider.setRange(8, 20)
        self.font_slider.setValue(current_size)
        
        self.theme_combo = QComboBox()
        self.theme_combo.addItems(["Dark", "Light"])
        self.theme_combo.setCurrentText(current_theme)
        
        form.addRow("Font Size:", self.font_slider)
        form.addRow("Theme:", self.theme_combo)
        layout.addLayout(form)
        
        btn_close = QPushButton("Apply & Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)

    def update_global_font(self):
        # Dark/Light Theme Mapping
        if self.theme == "Dark":
            bg, fg, btn, btn_hover = "#2d2d2d", "#ffffff", "#3e3e42", "#505050"
        else:
            bg, fg, btn, btn_hover = "#f0f0f0", "#000000", "#e0e0e0", "#d0d0d0"

        app = QApplication.instance()
        app.setStyleSheet(f"""
            QWidget {{ font-size: {self.current_font_size}pt; background-color: {bg}; color: {fg}; }}
            QPushButton {{ padding: 6px; background-color: {btn}; border: 1px solid #888; border-radius: 4px; }}
            QPushButton:hover {{ background-color: {btn_hover}; }}
            QLineEdit, QComboBox {{ background-color: #ffffff; color: #000000; }}
            QGroupBox {{ font-weight: bold; border: 1px solid #888; margin-top: 10px; padding-top: 10px; }}
        """)

# DATA MANAGERS (TABLE EDITORS)
class DataManagerDialog(QDialog):
    def __init__(self, parent, title, headers, data_list):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(600, 400)
        self.data_list = data_list
        self.headers = headers
        
        layout = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        btn_layout = QHBoxLayout()
        self.btn_add = QPushButton("+ Add Row")
        self.btn_remove = QPushButton("- Remove Selected")
        self.btn_apply = QPushButton("Apply & Close")
        
        self.btn_add.clicked.connect(self.add_row)
        self.btn_remove.clicked.connect(self.remove_row)
        self.btn_apply.clicked.connect(self.apply_and_close)
        
        btn_layout.addWidget(self.btn_add)
        btn_layout.addWidget(self.btn_remove)
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_apply)
        
        layout.addWidget(self.table)
        layout.addLayout(btn_layout)
        self.populate_table()

    def populate_table(self):
        self.table.setRowCount(len(self.data_list))
        for row_idx, row_data in enumerate(self.data_list):
            for col_idx, value in enumerate(row_data):
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))

    def add_row(self):
        row_position = self.table.rowCount()
        self.table.insertRow(row_position)
        for i in range(self.table.columnCount()):
            default_val = "Inlet" if self.headers[i] == "Type" else "0"
            self.table.setItem(row_position, i, QTableWidgetItem(default_val))

    def remove_row(self):
        current_row = self.table.currentRow()
        if current_row >= 0: self.table.removeRow(current_row)

    def apply_and_close(self):
        self.data_list.clear()
        for row in range(self.table.rowCount()):
            row_data = []
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                val = item.text() if item else ""
                if self.headers[col] != "Type":
                    try: val = float(val)
                    except ValueError: val = 0.0
                row_data.append(val)
            self.data_list.append(row_data)
        self.accept()

# MAIN APPLICATION ENGINE
class FVMSetupGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        # 1. Get screen resolution
        screen = QApplication.primaryScreen()
        screen_geometry = screen.availableGeometry()
        width = screen_geometry.width()
        height = screen_geometry.height()

        self.setWindowTitle("FVM Simulation IDE")
        # self.resize(1400, 900)
        self.resize(int(width * 0.8), int(height * 0.8))
        self.setMinimumSize(800, 600)

        
        # Initialize variables
        self.obstacles = []
        self.bcs = []
        self.current_font_size = 10
        self.theme = "Dark"
        self.cycles_data = []
        self.errors_data = []
        self.draw_mode_active = False 
        self.temp_points = []    

        # Background Process Management
        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.read_solver_output)
        self.process.finished.connect(self.solver_finished)

        # Initialize UI Components
        self.input_refine = QLineEdit("1") # Define this before adding to layout
        
        self.setup_ui()
        # self.canvas.scene().sigMouseClicked.connect(self.canvas_clicked)
        self.setup_menu_bar()
        self.update_global_font()
        self.master_render()
    
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # --- LEFT CONTROLS PANEL ---
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        control_container = QWidget()
        control_panel = QVBoxLayout(control_container)
        
        # 1. Domain Setup
        domain_group = QGroupBox("1. Domain Definition")
        domain_layout = QFormLayout()
        
        # Now these are purely data fields, not used for drawing
        self.input_lx = QLineEdit("3.0") 
        self.input_ly = QLineEdit("1.0") 
        
        # These will be used for drawing bounds
        self.input_nx = QLineEdit("120")    
        self.input_ny = QLineEdit("40")     
        self.input_refine = QLineEdit("1") # New refinement option
        self.btn_draw_bbox = QPushButton("Apply Grid Bounds")
        self.btn_draw_bbox.clicked.connect(self.master_render)
        # self.btn_draw_bbox.clicked.connect(self.auto_scale_view)
        
        domain_layout.addRow("Phys Length X (Lx):", self.input_lx)
        domain_layout.addRow("Phys Length Y (Ly):", self.input_ly)
        domain_layout.addRow("Grid Cells X (nx):", self.input_nx)
        domain_layout.addRow("Grid Cells Y (ny):", self.input_ny)
        domain_layout.addRow("No. of Refinements:", self.input_refine)
        domain_layout.addRow(self.btn_draw_bbox)
        domain_group.setLayout(domain_layout)

        # 2. Obstacles Definition
        obs_group = QGroupBox("2. Obstacles")
        obs_layout = QVBoxLayout()
        self.lbl_obs_status = QLabel("0 obstacles loaded.")
        
        offset_form = QFormLayout()
        self.input_obs_x_off = QLineEdit("0")
        self.input_obs_y_off = QLineEdit("0")
        self.input_obs_x_off.textChanged.connect(self.master_render)
        self.input_obs_y_off.textChanged.connect(self.master_render)
        offset_form.addRow("X Offset (Indices):", self.input_obs_x_off)
        offset_form.addRow("Y Offset (Indices):", self.input_obs_y_off)
        obs_layout.addLayout(offset_form)
        
        btn_layout_obs = QHBoxLayout()
        self.btn_browse_obs = QPushButton("Load .txt")
        self.btn_browse_obs.clicked.connect(self.load_obstacles_txt)
        self.btn_manage_obs = QPushButton("Manage List")
        self.btn_manage_obs.clicked.connect(self.open_obstacle_manager)
        btn_layout_obs.addWidget(self.btn_browse_obs)
        btn_layout_obs.addWidget(self.btn_manage_obs)
        obs_layout.addWidget(self.lbl_obs_status)
        obs_layout.addLayout(btn_layout_obs)
        obs_group.setLayout(obs_layout)

        # 3. Boundary Conditions
        bc_group = QGroupBox("3. Boundary Conditions")
        bc_layout = QVBoxLayout()
        self.lbl_bc_status = QLabel("0 BCs loaded.")
        btn_layout_bc = QHBoxLayout()
        self.btn_browse_bc = QPushButton("Load .txt")
        self.btn_browse_bc.clicked.connect(self.load_bc_txt)
        self.btn_manage_bc = QPushButton("Manage List")
        self.btn_manage_bc.clicked.connect(self.open_bc_manager)
        btn_layout_bc.addWidget(self.btn_browse_bc)
        btn_layout_bc.addWidget(self.btn_manage_bc)
        
        draw_layout = QHBoxLayout()
        self.bc_combo = QComboBox()
        self.bc_combo.addItems(["Inlet", "Outlet"])
        self.btn_draw_canvas = QPushButton("Draw Canvas Line")
        self.btn_draw_canvas.setCheckable(True)
        self.btn_draw_canvas.clicked.connect(self.toggle_draw_mode)
        draw_layout.addWidget(self.bc_combo)
        draw_layout.addWidget(self.btn_draw_canvas)

        bc_layout.addWidget(self.lbl_bc_status)
        bc_layout.addLayout(btn_layout_bc)
        bc_layout.addWidget(QLabel("Interactive Drawing (Indices):"))
        bc_layout.addLayout(draw_layout)
        bc_group.setLayout(bc_layout)

        # 4. Solver Settings
        solver_settings_group = QGroupBox("4. Solver Settings")
        solver_settings_layout = QFormLayout()
        
        self.solver_type_combo = QComboBox()
        self.solver_type_combo.addItems(["SIMPLE", "COUPLED"])
        # self.solver_type_combo.currentTextChanged.connect(self.toggle_multigrid)
        self.press_solver_combo = QComboBox()
        self.press_solver_combo.addItems(["Jacobi", "Gauss Seidel"])
        solver_settings_layout.addRow("Algorithm:", self.solver_type_combo)
        solver_settings_layout.addRow("Pressure Solver:", self.press_solver_combo)
        solver_settings_group.setLayout(solver_settings_layout)

        # 5. Execution Controls
        solver_group = QGroupBox("5. Execution Engine")
        solver_layout = QVBoxLayout()
        self.hardware_combo = QComboBox()
        self.hardware_combo.addItems(["CPU Execution (GCC)", "GPU Execution (OpenACC / NVC)"])
        self.btn_run_solver = QPushButton("Run Simulation Engine")
        self.btn_run_solver.setStyleSheet("background-color: #aaddaa; font-weight: bold;")
        self.btn_run_solver.clicked.connect(self.execute_solver_pipeline)
        
        solver_layout.addWidget(QLabel("Target Processing Unit:"))
        solver_layout.addWidget(self.hardware_combo)
        solver_layout.addWidget(self.btn_run_solver)
        solver_group.setLayout(solver_layout)

        # 5. Visualization Controls
        vis_group = QGroupBox("5. Visualization")
        vis_layout = QFormLayout()
        self.vis_mode_combo = QComboBox()
        self.vis_mode_combo.addItems(["None", "Streamlines", "Contours"])
        self.vis_mode_combo.currentTextChanged.connect(self.master_render)
        self.vis_var_combo = QComboBox()
        self.vis_var_combo.addItems(["ux", "uy", "umag", "p"])
        self.vis_var_combo.currentTextChanged.connect(self.master_render)
        vis_layout.addRow("Mode:", self.vis_mode_combo)
        vis_layout.addRow("Variable:", self.vis_var_combo)
        vis_group.setLayout(vis_layout)

        # ADD GROUPS TO PANEL
        control_panel.addWidget(domain_group)
        control_panel.addWidget(obs_group)
        control_panel.addWidget(bc_group)
        control_panel.addWidget(solver_settings_group) # Fixed: Now it is added to the layout
        control_panel.addWidget(solver_group)
        control_panel.addWidget(vis_group)
        control_panel.addStretch()
        scroll_area.setWidget(control_container)

        # 2. Prepare your Right Side (Canvas + Monitors)
        # Combine the canvas and bottom panels into one widget
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        self.canvas = pg.PlotWidget(title="Domain Geometry")
        self.canvas.setAspectLocked(True)

        # Bottom monitor area
        bottom_monitor = QWidget()
        bottom_layout = QHBoxLayout(bottom_monitor)
        self.terminal = QPlainTextEdit()
        self.convergence_plot = pg.PlotWidget(title="Residuals")
        self.convergence_curve = self.convergence_plot.plot(pen=pg.mkPen('g', width=2))
        self.convergence_curve.setData([], [])
        bottom_layout.addWidget(self.terminal)
        bottom_layout.addWidget(self.convergence_plot)
        
        right_layout.addWidget(self.canvas)
        right_layout.addWidget(bottom_monitor)
        QTimer.singleShot(100, lambda: self.canvas.scene().sigMouseClicked.connect(self.canvas_clicked))

        # --- SPLITTER LAYOUT ---
        main_splitter = QSplitter(Qt.Horizontal)
        main_splitter.addWidget(scroll_area)
        main_splitter.addWidget(right_panel)
        # Set splitter to take 20% for the left panel and 80% for the right
        splitter_width = int(self.width())
        main_splitter.setSizes([int(splitter_width * 0.25), int(splitter_width * 0.75)])

        # 4. Add the splitter to the main layout
        main_layout.addWidget(main_splitter)

    # --- MENU CONFIGURATION ---
    def setup_menu_bar(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        edit_menu = menubar.addMenu("Edit")
        
        exit_action = file_menu.addAction("Exit")
        exit_action.triggered.connect(self.close)
        
        pref_action = edit_menu.addAction("Preferences")
        pref_action.triggered.connect(self.open_preferences)

    def open_preferences(self):
        dialog = PreferencesDialog(self, self.current_font_size, self.theme)
        if dialog.exec():
            self.current_font_size = dialog.font_slider.value()
            self.theme = dialog.theme_combo.currentText()
            self.update_global_font()

    def update_global_font(self):
        # Use the class attribute self.current_font_size
        size = self.current_font_size
        
        # Theme logic
        if self.theme == "Dark":
            bg, fg, btn, btn_hover = "#2d2d2d", "#ffffff", "#3e3e42", "#505050"
        else:
            bg, fg, btn, btn_hover = "#f0f0f0", "#000000", "#e0e0e0", "#d0d0d0"

        app = QApplication.instance()
        app.setStyleSheet(f"""
            QWidget {{ font-size: {size}pt; background-color: {bg}; color: {fg}; }}
            QGroupBox {{ font-weight: bold; border: 1px solid #888; margin-top: 10px; padding-top: 10px; }}
            QPushButton {{ padding: 6px; background-color: {btn}; border: 1px solid #888; border-radius: 4px; color: {fg}; }}
            QPushButton:hover {{ background-color: {btn_hover}; }}
            QLineEdit, QComboBox {{ background-color: #ffffff; color: #000000; }}
        """)

    # --- DIALOG CALLS ---
    def open_obstacle_manager(self):
        headers = ["i_init", "i_final", "j_init", "j_final"]
        if DataManagerDialog(self, "Manage Obstacles", headers, self.obstacles).exec():
            self.lbl_obs_status.setText(f"{len(self.obstacles)} obstacles loaded.")
            self.master_render()

    def open_bc_manager(self):
        headers = ["Type", "i_init", "i_final", "j_init", "j_final", "u", "v"]
        if DataManagerDialog(self, "Manage Boundaries", headers, self.bcs).exec():
            self.lbl_bc_status.setText(f"{len(self.bcs)} BCs loaded.")
            self.master_render()

    # --- INPUT PROCESSING PARSERS ---
    def load_obstacles_txt(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Obstacle File", "", "Text Files (*.txt)")
        if not file_path: return
        try:
            with open(file_path, 'r') as file:
                lines = [line.strip() for line in file.readlines() if line.strip()]
            self.obstacles.clear()
            for i in range(1, len(lines)):
                parts = lines[i].replace(',', ' ').split()
                if len(parts) >= 4:
                    self.obstacles.append([float(p) for p in parts[:4]])
            self.lbl_obs_status.setText(f"{len(self.obstacles)} obstacles loaded.")
            self.master_render()
        except Exception as e: QMessageBox.critical(self, "Error", f"Failed to parse.\n{e}")

    def load_bc_txt(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select BC File", "", "Text Files (*.txt)")
        if not file_path: return
        try:
            with open(file_path, 'r') as file:
                lines = [line.strip() for line in file.readlines() if line.strip()]
            self.bcs.clear()
            idx = 0
            n_inlets = int(lines[idx]); idx += 1
            for _ in range(n_inlets):
                parts = lines[idx].replace(',', ' ').split()
                if len(parts) >= 4:
                    self.bcs.append(["Inlet", float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), 0.0, 0.0])
                idx += 1
            n_outlets = int(lines[idx]); idx += 1
            for _ in range(n_outlets):
                parts = lines[idx].replace(',', ' ').split()
                if len(parts) >= 4:
                    self.bcs.append(["Outlet", float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), 0.0, 0.0])
                idx += 1
            self.lbl_bc_status.setText(f"{len(self.bcs)} BCs loaded.")
            self.master_render()
        except Exception as e: QMessageBox.critical(self, "Error", f"Failed to parse BC file.\n{e}")


    # --- SOLVER EXECUTION PIPELINE ---
    def execute_solver_pipeline(self):
        nx = self.input_nx.text()
        ny = self.input_ny.text()
        ngrid = self.input_refine.text()
        lx = self.input_lx.text()
        ly = self.input_ly.text()
        
        try:
            x_off = float(self.input_obs_x_off.text())
            y_off = float(self.input_obs_y_off.text())
        except ValueError:
            x_off, y_off = 0.0, 0.0

        try:
            # pattern.txt
            with open("pattern.txt", "w") as f:
                f.write(f"{len(self.obstacles)}\n")
                for obs in self.obstacles:
                    f.write(f"{int(obs[0])} {int(obs[1])} {int(obs[2])} {int(obs[3])}\n")

            # input_data.txt
            with open("inputData2D", "w") as f:
                f.write(f"{nx} {ny} {ngrid} 200 25 25 2 1.0 .00001\n")
                f.write(f"{lx} {ly} .01 1.0 0.15 0.15\n")
                f.write("0 0.0 0 0.0\n")
                f.write(f"{len(self.obstacles)} {x_off} {y_off}\n")
                
                inlets = [b for b in self.bcs if b[0] == "Inlet"]
                outlets = [b for b in self.bcs if b[0] == "Outlet"]
                f.write(f"{len(inlets)} {len(outlets)}\n")
                
                for b in inlets:
                    f.write(f"{b[3]} {b[4]} {b[5]} {b[6]}\n")
                for b in outlets:
                    f.write(f"{b[3]} {b[4]} {b[5]} {b[6]}\n")
                    
        except Exception as e:
            QMessageBox.critical(self, "File System Error", f"Failed to write parameter assets:\n{e}")
            return

        self.terminal.clear()
        self.cycles_data.clear()
        self.errors_data.clear()
        self.convergence_curve.setData([], [])

        is_gpu = "GPU" in self.hardware_combo.currentText()
        if is_gpu:
            compile_cmd = "nvc -acc -Minfo=accel FVMsolver.c -o FVMsolver.exe"
            run_cmd = "FVMsolver.exe"
        else:
            compile_cmd = "gcc -O3 FVMsolver.c -o FVMsolver.exe"
            run_cmd = "FVMsolver.exe"

        self.terminal.appendPlainText(f"Executing Configuration: {compile_cmd}\n")
        
        ret = os.system(compile_cmd)
        if ret != 0:
            self.terminal.appendPlainText("[CRITICAL] Compilation sequence failure.\n")
            return
            
        self.terminal.appendPlainText("[SUCCESS] Structural binaries built. Running Model...\n")
        self.process.start(run_cmd)

    def render_vtk_data(self, mode):
        if not os.path.exists("vtk_acc.vtk"): return
        
        try:
            mesh = pv.read("vtk_acc.vtk")
            refine = int(self.input_refine.text())
            scale = 1.0 * float(self.input_nx.text()) / (2**(refine-1) * float(self.input_lx.text()))

            # --- DYNAMIC DATA DERIVATION ---
            # Extract Velocity vector components (PyVista stores them as 'Velocity')
            # The 'Velocity' array is a (N, 3) array.
            vel = mesh.point_data["Velocity"] 
            mesh["ux"] = vel[:, 0]
            mesh["uy"] = vel[:, 1]
            # Calculate magnitude on the fly
            mesh["umag"] = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2)

            if mode == "Streamlines":
                streams = mesh.streamlines(vectors="Velocity", n_points=30)
                for i in range(streams.n_lines):
                    line = streams.get_cell(i)
                    pts = line.points * scale
                    self.canvas.plot(pts[:, 0], pts[:, 1], pen=pg.mkPen('y', width=1))

            elif mode == "Contours":
                var = self.vis_var_combo.currentText()
                # 'P' is already in mesh.array_names
                if var in mesh.array_names or var in ["ux", "uy", "umag"]:
                    contours = mesh.contour(isosurfaces=8, scalars=var)
                    for i in range(contours.n_cells):
                        cell = contours.get_cell(i)
                        pts = cell.points * scale
                        self.canvas.plot(pts[:, 0], pts[:, 1], pen=pg.mkPen('c', width=1))
                    
        except Exception as e:
            self.terminal.appendPlainText(f"[ERROR] Visualization failed: {e}")
    
    def read_solver_output(self):
        data = self.process.readAllStandardOutput().data().decode()
        self.terminal.appendPlainText(data)
        
        pattern = r"icycle\s*=\s*(\d+)\s+error\s*=\s*([\d\.eE\+\-]+)"
        for line in data.splitlines():
            match = re.search(pattern, line)
            if match:
                self.cycles_data.append(int(match.group(1)))
                self.errors_data.append(float(match.group(2)))
                
        if self.cycles_data:
            self.convergence_curve.setData(self.cycles_data, self.errors_data)

    def solver_finished(self, exit_code, exit_status):
        self.terminal.appendPlainText(f"\n[TERMINATED] Status code: {exit_code}")

    # --- INTERACTIVE CANVAS LOGIC ---
    def toggle_draw_mode(self, checked):
        self.draw_mode_active = checked
        self.temp_points = []
        if checked:
            self.btn_draw_canvas.setText("Click 2 Viewport Pts")
            self.canvas.setCursor(Qt.CrossCursor)
            self.canvas.setMouseEnabled(x=False, y=False)
        else:
            self.btn_draw_canvas.setText("Draw Canvas Line")
            self.canvas.setCursor(Qt.ArrowCursor)
            self.canvas.setMouseEnabled(x=True, y=True)
            self.master_render()

    def canvas_clicked(self, event):
        if not hasattr(self, 'draw_mode_active') or not self.draw_mode_active:
            return
        if not self.draw_mode_active or event.button() != Qt.LeftButton: return
        pos = event.scenePos()
        view_pos = self.canvas.getPlotItem().vb.mapSceneToView(pos)
        snap_x, snap_y = round(view_pos.x()), round(view_pos.y())
        self.temp_points.append((snap_x, snap_y))
        
        self.canvas.addItem(pg.ScatterPlotItem([snap_x], [snap_y], size=10, pen='y', brush='y'))
        
        if len(self.temp_points) == 2:
            bc_type = self.bc_combo.currentText()
            xi, yi = self.temp_points[0]
            xf, yf = self.temp_points[1]
            self.bcs.append([bc_type, xi, xf, yi, yf, 0.0, 0.0])
            self.lbl_bc_status.setText(f"{len(self.bcs)} BCs loaded.")
            self.temp_points = []
            self.master_render()

    def auto_scale_view(self):
        try:
            nx = float(self.input_nx.text())
            ny = float(self.input_ny.text())
            self.canvas.setXRange(-nx*0.1, nx*1.1, padding=0)
            self.canvas.setYRange(-ny*0.1, ny*1.1, padding=0)
        except ValueError: pass

    def master_render(self):
        """Single rendering loop based purely on integer indices."""
        self.canvas.clear()
        self.canvas.showGrid(x=True, y=True)
        self.canvas.setAspectLocked(True)

        # 1. Base Domain uses grid counts (nx, ny) instead of physical lengths
        try:
            nx = float(self.input_nx.text())
            ny = float(self.input_ny.text())
            bbox = QGraphicsRectItem(0, 0, nx, ny)
            bbox.setPen(pg.mkPen('w', width=3)) 
            self.canvas.addItem(bbox)
        except ValueError: pass

        # 2. Structural Offset Evaluation
        try:
            x_off = float(self.input_obs_x_off.text())
            y_off = float(self.input_obs_y_off.text())
        except ValueError: 
            x_off, y_off = 0.0, 0.0

        # 3. Draw Obstacles directly as indices
        for obs in self.obstacles:
            i_init, i_final, j_init, j_final = obs
            
            xi = i_init + x_off
            xf = i_final + x_off
            yi = j_init + y_off
            yf = j_final + y_off
            
            rect = QGraphicsRectItem(xi, yi, xf - xi, yf - yi)
            rect.setPen(pg.mkPen((150, 150, 150), width=1)) 
            rect.setBrush(pg.mkBrush(80, 80, 80, 200)) 
            self.canvas.addItem(rect)

        # 4. Draw Boundaries
        for bc in self.bcs:
            bc_type, xi, xf, yi, yf, u, v = bc
            color = (50, 150, 255) if "Inlet" in bc_type else (255, 50, 50)
            self.canvas.addItem(pg.PlotCurveItem([xi, xf], [yi, yf], pen=pg.mkPen(color, width=5)))

        # 5. Render Results if available
        mode = self.vis_mode_combo.currentText()
        if mode != "None":
            self.render_vtk_data(mode)


if __name__ == "__main__":
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QApplication(sys.argv)
    window = FVMSetupGUI()
    window.show()
    sys.exit(app.exec())