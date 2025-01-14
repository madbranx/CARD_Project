class MaterialProperty:
    def __init__(self, collision_diameter):
        self.collision_diameter = collision_diameter

    def density(self, temperature):
        # Placeholder function for density as a function of temperature
        # Replace with actual formula or data
        return 1000 - 0.1 * (temperature - 273.15)  # Example: linear decrease with temperature

    def heat_capacity(self, temperature):
        # Placeholder function for heat capacity as a function of temperature
        # Replace with actual formula or data
        return 4.18  # Example: constant value in J/(g*K)

    def heat_conductivity(self, temperature):
        # Placeholder function for heat conductivity as a function of temperature
        # Replace with actual formula or data
        return 0.6 + 0.001 * (temperature - 273.15)  # Example: linear increase with temperature
