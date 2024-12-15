import tkinter as tk
from tkinter import simpledialog

class DisjointSet:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)
        if rootX != rootY:
            if self.rank[rootX] > self.rank[rootY]:
                self.parent[rootY] = rootX
            elif self.rank[rootX] < self.rank[rootY]:
                self.parent[rootX] = rootY
            else:
                self.parent[rootY] = rootX
                self.rank[rootX] += 1

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0


def rename_columns_with_dsu(first_array, second_array):
    dsu = DisjointSet()

    for pair in second_array:
        for element in pair:
            dsu.add(element)

    for pair in second_array:
        dsu.union(pair[0], pair[1])
    
    root_ground = dsu.find('Ground')
    dsu.parent['Ground'] = 'Ground'
    dsu.parent[root_ground] = 'Ground'
    
    root_to_number = {}
    current_number = 0
    for element in dsu.parent:
        root = dsu.find(element)
        if root not in root_to_number:
            if root == "Ground":
                root_to_number[root] = 0 
            else:
                current_number += 1
                root_to_number[root] = current_number

    renamed_array = []
    for row in first_array:
        renamed_row = row[:]
        renamed_row[1] = str(root_to_number[dsu.find(row[1])]) 
        renamed_row[2] = str(root_to_number[dsu.find(row[2])]) 
        renamed_array.append(renamed_row)

    return renamed_array

class Component:
    def __init__(self, component_id, terminals, value=None):
        self.component_id = component_id
        self.terminals = terminals
        self.value = value

class CircuitGraph:
    def __init__(self):
        self.components = []
        self.connections = [] 

    def add_component(self, component):
        self.components.append(component)

    def add_connection(self, comp1_id, term1, comp2_id, term2):
        self.connections.append((comp1_id, term1, comp2_id, term2))

    def generate_netlist(self):
        netlist = []
        nodal=[]
        final_netlist=[]
        final_nodal=[]
        ff_netlist=[]
        for component in self.components:
            if (component.component_id[0]!="G"):
                entry = f"{component.component_id} " + " ".join(component.terminals)
                if component.value:
                    entry += f" {component.value}"
                netlist.append(entry)
        
        for comp1, term1, comp2, term2 in self.connections:
            nodal.append(f"{term1} connected to {term2}")
        
        for entry in netlist:
            final_netlist.append(entry.split(" "))
            
        for entry in nodal:
            final_nodal.append(entry.split(" connected to "))
            
        ff_netlist=rename_columns_with_dsu(final_netlist, final_nodal)
        return ff_netlist

class CircuitGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GSpice V_2.0.4")
        self.iconbitmap("logo.ico")
        self.geometry("800x600")
        
        self.canvas = tk.Canvas(self, bg="white")
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.canvas.create_text(65,10,text="Wiring option is OFF",tags="off_condition")

        self.add_resistor_button = tk.Button(self, text="Resistor", command=self.select_resistor)
        self.add_resistor_button.pack(side=tk.LEFT)

        self.add_capacitor_button = tk.Button(self, text="Capacitor", command=self.select_capacitor)
        self.add_capacitor_button.pack(side=tk.LEFT)

        self.add_inductor_button = tk.Button(self, text="Inductor", command=self.select_inductor)
        self.add_inductor_button.pack(side=tk.LEFT)

        self.add_voltage_source_button = tk.Button(self, text="Voltage Source", command=self.select_voltage_source)
        self.add_voltage_source_button.pack(side=tk.LEFT)
        
        self.add_current_source_button = tk.Button(self, text="Current Source", command=self.select_current_source)
        self.add_current_source_button.pack(side=tk.LEFT)

        self.add_ground_button = tk.Button(self, text="Ground", command=self.select_ground)
        self.add_ground_button.pack(side=tk.LEFT)

        self.add_wire_button = tk.Button(self, text="Wire", command=self.toggle_wire_mode)
        self.add_wire_button.pack(side=tk.LEFT)

        self.generate_simulate_button = tk.Button(self, text="Simulate", command=self.simulate)
        self.generate_simulate_button.pack(side=tk.RIGHT)

        self.circuit = CircuitGraph()
        self.component_counter = 1
        self.selected_component_type = None
        self.wire_mode = False
        self.nodes = {} 
        self.selected_terminals = [] 
        
        self.canvas.bind("<Button-1>", self.place_component)

    def select_resistor(self):
        self.selected_component_type = "resistor"
        self.wire_mode = False

    def select_capacitor(self):
        self.selected_component_type = "capacitor"
        self.wire_mode = False

    def select_inductor(self):
        self.selected_component_type = "inductor"
        self.wire_mode = False

    def select_voltage_source(self):
        self.selected_component_type = "voltage_source"
        self.wire_mode = False
        
    def select_current_source(self):
        self.selected_component_type = "current_source"
        self.wire_mode = False

    def select_ground(self):
        self.selected_component_type = "ground"
        self.wire_mode = False

    def toggle_wire_mode(self):
        self.wire_mode = not self.wire_mode
        self.selected_component_type = "wire" if self.wire_mode else None
        if (self.wire_mode):
            self.canvas.delete("node_condition")
            self.canvas.delete("off_condition")
            self.canvas.create_text(65,10,text="Wiring option is ON",tags="on_condition")
        else:
            self.canvas.delete("node_condition")
            self.canvas.delete("on_condition")
            self.canvas.create_text(65,10,text="Wiring option is OFF",tags="off_condition")            

    def place_component(self, event):
        if not self.selected_component_type or self.selected_component_type == "wire":
            return

        if (self.selected_component_type=="inductor"):
            component_id = f"L{self.component_counter}"
            self.component_counter += 1
            value = None
        elif (self.selected_component_type=="current_source"):
            component_id = f"I{self.component_counter}"
            self.component_counter += 1
            value = None
        else:                
            component_id = f"{self.selected_component_type[0].upper()}{self.component_counter}"
            self.component_counter += 1
            value = None

        x, y = event.x, event.y

        if self.selected_component_type == "resistor":
            value = simpledialog.askstring("Resistor Value", "Enter resistor value (e.g., 1E3 for 1kohm):")
            if value is None:
                return 
            else:
                self.canvas.create_rectangle(x, y, x + 60, y + 20, fill="yellow", tags=component_id)
                self.canvas.create_text(x + 30, y + 10, text=component_id)
                terminals = {f"{component_id}.n1": (x, y + 10), f"{component_id}.n2": (x + 60, y + 10)}

        elif self.selected_component_type == "capacitor":
            value = simpledialog.askstring("Capacitor Value", "Enter capacitor value (e.g., 10E-9 for 10nF):")
            if value is None:
                return 
            else:
                self.canvas.create_rectangle(x, y, x + 60, y + 20, fill="cyan", tags=component_id)
                self.canvas.create_text(x + 30, y + 10, text=component_id)
                terminals = {f"{component_id}.n1": (x, y + 10), f"{component_id}.n2": (x + 60, y + 10)}

        elif self.selected_component_type == "inductor":
            value = simpledialog.askstring("Inductor Value", "Enter inductor value (e.g., 1E-3 for 1mH):")
            if value is None:
                return 
            else:
                self.canvas.create_rectangle(x, y, x + 60, y + 20, fill="orange", tags=component_id)
                self.canvas.create_text(x + 30, y + 10, text=component_id)
                terminals = {f"{component_id}.n1": (x, y + 10), f"{component_id}.n2": (x + 60, y + 10)}

        elif self.selected_component_type == "voltage_source":
            value = simpledialog.askstring("Voltage Source Value", "Enter voltage value (e.g., 5 for 5V):")
            if value is None:
                return 
            else:
                self.canvas.create_rectangle(x, y, x + 60, y + 20, fill="green", tags=component_id)
                self.canvas.create_text(x + 30, y + 10, text=component_id)
                terminals = {f"{component_id}.n1": (x, y + 10), f"{component_id}.n2": (x + 60, y + 10)}
                
        elif self.selected_component_type == "current_source":
            value = simpledialog.askstring("Current Source Value", "Enter current value (e.g., 5 for 5A):")
            if value is None:
                return 
            else:
                self.canvas.create_rectangle(x, y, x + 60, y + 20, fill="red", tags=component_id)
                self.canvas.create_text(x + 30, y + 10, text=component_id)
                terminals = {f"{component_id}.n1": (x, y + 10), f"{component_id}.n2": (x + 60, y + 10)}

        elif self.selected_component_type == "ground":
            self.canvas.create_oval(x - 10, y - 10, x + 10, y + 10, fill="brown", tags=component_id)
            self.canvas.create_text(x, y + 20, text=component_id)
            terminals = {"Ground": (x, y)}

        component = Component(component_id, terminals, value)
        self.circuit.add_component(component)
        self.nodes[component_id] = terminals 
       
        for term, (tx, ty) in terminals.items():
            terminal_tag = f"{component_id}_{term}"
            self.canvas.create_oval(tx - 3, ty - 3, tx + 3, ty + 3, fill="black", tags=terminal_tag)
            
            self.canvas.tag_bind(terminal_tag, "<Button-1>", lambda event, cid=component_id, t=term: self.select_terminal(cid, t))

        self.selected_component_type = None

    def select_terminal(self, component_id, terminal):
        if self.wire_mode:
            if len(self.selected_terminals)==0:
                self.canvas.delete("node_condition")
                self.canvas.create_text(65,25,text="Select first node",tags="node_condition")
                
            if (component_id, terminal) not in self.selected_terminals:
                self.selected_terminals.append((component_id, terminal))
            
            if len(self.selected_terminals)==1:
                self.canvas.delete("node_condition")
                self.canvas.create_text(65,25,text="Select second node",tags="node_condition")
            
            if len(self.selected_terminals) == 2:
                (comp1, term1), (comp2, term2) = self.selected_terminals
                x1, y1 = self.nodes[comp1][term1]
                x2, y2 = self.nodes[comp2][term2]

                self.canvas.create_line(x1, y1, x2, y2, fill="black", width=2)

                self.circuit.add_connection(comp1, term1, comp2, term2)

                self.selected_terminals = []
        
    def simulate(self):
        netlist = self.circuit.generate_netlist()
        with open("output.txt", "w") as file:
            for entry in netlist:
                for component in entry:
                    print(component,end=" ")
                    file.write(str(component) + " ")
                print("\n")
                file.write("\n")
        file.close()
        
if __name__ == "__main__":
    app = CircuitGUI()
    app.mainloop()