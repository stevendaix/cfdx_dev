"""Contract-driven Qt setup editor for CFDX cases."""
from __future__ import annotations
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox,QComboBox,QDoubleSpinBox,QFormLayout,QHBoxLayout,QLineEdit,QListWidget,QPushButton,QTabWidget,QVBoxLayout,QWidget
from .boundary_setup import SCALAR_TYPES, VELOCITY_TYPES, FIELD_SPECS
from .case import Case
from .initialization import InitializationMode, InitializationSpec
from .materials import MaterialSpec
from .physics_setup import PHYSICS_SPECS,TURBULENCE_MODELS

class CaseSetupPanel(QWidget):
    changed=Signal()
    """Thin Qt adapter over the shared headless setup contracts."""

    def __init__(self,case: Case,parent: QWidget|None=None)->None:
        super().__init__(parent); self.case=case
        tabs=QTabWidget(self)
        tabs.addTab(self._physics_tab(),"Physics"); tabs.addTab(self._material_tab(),"Materials")
        tabs.addTab(self._boundary_tab(),"Boundaries"); tabs.addTab(self._initialization_tab(),"Initialization")
        layout=QVBoxLayout(self); layout.addWidget(tabs)

    @staticmethod
    def _real(value: float)->QDoubleSpinBox:
        box=QDoubleSpinBox(); box.setRange(-1e300,1e300); box.setDecimals(10); box.setValue(value); return box

    def _physics_tab(self)->QWidget:
        w=QWidget(); form=QFormLayout(w)
        self.physics_model=QComboBox(); self.physics_model.addItems([s.key for s in PHYSICS_SPECS])
        self.physics_enabled=QCheckBox("Enabled"); self.turbulence_model=QComboBox(); self.turbulence_model.addItems(TURBULENCE_MODELS)
        self.physics_button=QPushButton("Apply physics")
        form.addRow("Model",self.physics_model); form.addRow("",self.physics_enabled); form.addRow("Turbulence model",self.turbulence_model); form.addRow("",self.physics_button)
        self.physics_model.currentTextChanged.connect(lambda m:self.turbulence_model.setEnabled(m=="turbulence"))
        self.physics_button.clicked.connect(self._apply_physics); self.turbulence_model.setEnabled(False)
        return w

    def _material_tab(self)->QWidget:
        w=QWidget(); form=QFormLayout(w)
        self.material_name=QLineEdit(); self.material_density=self._real(1.0); self.material_viscosity=self._real(1e-3)
        self.material_cp=self._real(1000.0); self.material_conductivity=self._real(0.0); self.material_button=QPushButton("Apply material")
        form.addRow("Name",self.material_name); form.addRow("Density [kg/m³]",self.material_density); form.addRow("Dynamic viscosity [Pa·s]",self.material_viscosity)
        form.addRow("Cp [J/(kg·K)]",self.material_cp); form.addRow("Conductivity [W/(m·K)]",self.material_conductivity); form.addRow("",self.material_button)
        self.material_button.clicked.connect(self._apply_material); return w

    def _boundary_tab(self)->QWidget:
        w=QWidget(); form=QFormLayout(w)
        self.boundary_list=QListWidget(); self.boundary_name=QLineEdit(); self.boundary_type=QComboBox()
        self.boundary_type.addItems(["inlet","outlet","wall","symmetry","periodic","interface","empty"])
        self.boundary_field=QComboBox(); self.boundary_field.addItems(list(FIELD_SPECS))
        self.boundary_scalar_type=QComboBox(); self.boundary_scalar_type.addItems(SCALAR_TYPES)
        self.boundary_value=self._real(0.0); self.boundary_gradient=self._real(0.0)
        self.boundary_x=self._real(0.0); self.boundary_y=self._real(0.0); self.boundary_z=self._real(0.0)
        self.boundary_button=QPushButton("Apply / Update"); self.boundary_remove_button=QPushButton("Remove selected")
        form.addRow("Existing patches",self.boundary_list); form.addRow("Name",self.boundary_name); form.addRow("Patch type",self.boundary_type)
        form.addRow("Field",self.boundary_field); form.addRow("Scalar condition",self.boundary_scalar_type); form.addRow("Value",self.boundary_value); form.addRow("Gradient",self.boundary_gradient)
        form.addRow("Vector X [m/s]",self.boundary_x); form.addRow("Vector Y [m/s]",self.boundary_y); form.addRow("Vector Z [m/s]",self.boundary_z)
        row=QHBoxLayout(); row.addWidget(self.boundary_button); row.addWidget(self.boundary_remove_button); form.addRow("",row)
        self.boundary_button.clicked.connect(self._apply_boundary); self.boundary_remove_button.clicked.connect(self._remove_boundary)
        self.boundary_list.currentTextChanged.connect(self._load_boundary); self.boundary_field.currentTextChanged.connect(self._update_boundary_field_controls); self._update_boundary_field_controls(self.boundary_field.currentText()); self._refresh_boundaries(); return w

    def _initialization_tab(self)->QWidget:
        w=QWidget(); form=QFormLayout(w); self.initialization_mode=QComboBox()
        self.initialization_mode.addItems([m.value for m in InitializationMode]); self.initialization_field=QLineEdit("U")
        self.initialization_value=self._real(0.0); self.initialization_button=QPushButton("Apply initialization")
        form.addRow("Mode",self.initialization_mode); form.addRow("Field",self.initialization_field); form.addRow("Uniform value",self.initialization_value); form.addRow("",self.initialization_button)
        self.initialization_button.clicked.connect(self._apply_initialization)
        self.initialization_mode.currentTextChanged.connect(lambda m:self.initialization_value.setEnabled(m=="uniform")); self.initialization_value.setEnabled(True)
        return w

    def set_case(self,case: Case)->None:
        self.case=case; self._refresh_boundaries()

    def _apply_physics(self)->None:
        model=self.physics_model.currentText(); values={"enabled":self.physics_enabled.isChecked()}
        if model=="turbulence": values["model"]=self.turbulence_model.currentText()
        self.case.physics[model]=values; self.changed.emit()

    def _apply_material(self)->None:
        name=self.material_name.text().strip()
        if not name:return
        material=MaterialSpec(name,self.material_density.value(),self.material_viscosity.value(),self.material_cp.value(),self.material_conductivity.value())
        try: material.validate()
        except ValueError:return
        self.case.materials[name]={"density":material.density,"dynamic_viscosity":material.dynamic_viscosity,"cp":material.cp,"conductivity":material.conductivity}; self.changed.emit()

    def _update_boundary_field_controls(self, field: str) -> None:
        vector = field == "velocity"
        for box in (self.boundary_x, self.boundary_y, self.boundary_z): box.setVisible(vector)
        self.boundary_value.setVisible(not vector); self.boundary_gradient.setVisible(not vector)
        self.boundary_scalar_type.clear(); self.boundary_scalar_type.addItems(VELOCITY_TYPES if vector else SCALAR_TYPES)

    def _apply_boundary(self)->None:
        name=self.boundary_name.text().strip()
        if not name:return
        values={"type":self.boundary_type.currentText()}
        if values["type"] not in {"periodic","interface","empty"}:
            field=self.boundary_field.currentText(); values["fields"]=[field]
            if field=="velocity": values.update({"velocity_type":self.boundary_scalar_type.currentText(),"velocity_value":[self.boundary_x.value(),self.boundary_y.value(),self.boundary_z.value()]})
            else: values.update({"scalar_type":self.boundary_scalar_type.currentText(),"value":self.boundary_value.value(),"gradient":self.boundary_gradient.value()})
        self.case.set_boundary(name,**values); self._refresh_boundaries(); self.changed.emit()

    def _remove_boundary(self)->None:
        item=self.boundary_list.currentItem()
        if item is not None: self.case.boundaries.pop(item.text(),None); self._refresh_boundaries(); self.changed.emit()

    def _refresh_boundaries(self)->None:
        self.boundary_list.blockSignals(True); self.boundary_list.clear(); self.boundary_list.addItems(list(self.case.boundaries)); self.boundary_list.blockSignals(False)

    def _load_boundary(self,name: str)->None:
        if not name:return
        b=self.case.boundaries.get(name,{})
        i=self.boundary_type.findText(str(b.get("type","wall"))); self.boundary_type.setCurrentIndex(max(0,i))
        i=self.boundary_scalar_type.findText(str(b.get("scalar_type","ZERO_GRADIENT"))); self.boundary_scalar_type.setCurrentIndex(max(0,i))
        field=str((b.get("fields") or ["pressure"])[0]); i=self.boundary_field.findText(field); self.boundary_field.setCurrentIndex(max(0,i))
        vector=b.get("velocity_value",[0.0,0.0,0.0]); self.boundary_x.setValue(float(vector[0])); self.boundary_y.setValue(float(vector[1])); self.boundary_z.setValue(float(vector[2]))
        self.boundary_value.setValue(float(b.get("value",0.0))); self.boundary_gradient.setValue(float(b.get("gradient",0.0))); self.boundary_name.setText(name)

    def _apply_initialization(self)->None:
        mode=InitializationMode(self.initialization_mode.currentText()); field=self.initialization_field.text().strip()
        spec=InitializationSpec(mode,field,self.initialization_value.value() if mode is InitializationMode.UNIFORM else None)
        try: spec.validate()
        except ValueError:return
        self.case.physics["initialization"]={"mode":mode.value,"field":spec.field,"value":spec.value}; self.changed.emit()
