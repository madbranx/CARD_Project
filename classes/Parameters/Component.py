import casadi as CasADi

class MaterialProperty:
    CONSTANT = 1
    LINEAR = 2
    POLYNOMIAL = 3

    def __init__(self, value=None, dependency_type=CONSTANT):
        self.value = value
        self.dependency_type = dependency_type

    def get_value(self, temperature=None):
        if self.dependency_type == self.CONSTANT:
            expr = self.value
        elif self.dependency_type == self.LINEAR:
            expr = self.value[0] * temperature + self.value[1]
        elif self.dependency_type == self.POLYNOMIAL:
            expr = sum(coefficient * (temperature ** i) for i, coefficient in enumerate(self.value))
        else:
            return None
        return expr

class Component:
    DENSITY = 1
    HEAT_CAPACITY = 2
    THERMAL_CONDUCTIVITY = 3
    COLLISION_AREA = 4
    DIFFUSION_VOLUME = 5
    DYNAMIC_VISCOSITY = 6
    MOLECULAR_WEIGHT = 7

    def __init__(self, name):
        self.name = name
        self.collision_area = None
        self.diffusion_volume = None
        self.density = None
        self.heat_capacity = None
        self.thermal_conductivity = None
        self.dynamic_viscosity = None
        self.molecular_weight = None

    def get_Name(self):
        return self.name

    def add_property(self, property_type, value=None, dependency_type=MaterialProperty.CONSTANT):
        if property_type == Component.DENSITY:
            self.density = MaterialProperty(value, dependency_type)

        elif property_type == Component.HEAT_CAPACITY:
            self.heat_capacity = MaterialProperty(value, dependency_type)

        elif property_type == Component.THERMAL_CONDUCTIVITY:
            self.thermal_conductivity = MaterialProperty(value, dependency_type)

        elif property_type == Component.COLLISION_AREA:
            self.collision_area = MaterialProperty(value, dependency_type)

        elif property_type == Component.DIFFUSION_VOLUME:
            self.diffusion_volume = MaterialProperty(value, dependency_type)

        elif property_type == Component.DYNAMIC_VISCOSITY:
            self.dynamic_viscosity = MaterialProperty(value, dependency_type)

        elif property_type == Component.MOLECULAR_WEIGHT:
            self.molecular_weight = MaterialProperty(value, dependency_type)


    def get_property(self, property_type, temperature=None):
        if property_type == Component.DENSITY:
            return self.density.get_value(temperature)
        elif property_type == Component.HEAT_CAPACITY:
            return self.heat_capacity.get_value(temperature)
        elif property_type == Component.THERMAL_CONDUCTIVITY:
            return self.thermal_conductivity.get_value(temperature)
        elif property_type == Component.COLLISION_AREA:
            return self.collision_area.get_value()
        elif property_type == Component.DIFFUSION_VOLUME:
            return self.diffusion_volume.get_value()
        elif property_type == Component.DYNAMIC_VISCOSITY:
            return self.dynamic_viscosity.get_value(temperature)
        elif property_type == Component.MOLECULAR_WEIGHT:
            return self.molecular_weight.get_value()

