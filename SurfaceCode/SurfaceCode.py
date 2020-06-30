# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 03:04:39 2020

@author: Andy
"""

from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit

class SurfaceCode():
    """
    Implementation of a distance d repetition code, implemented over
    T syndrome measurement rounds.
    """

    def __init__(self, d, T, mode='normal'):
        """
        Creates the circuits corresponding to a logical 0 and 1 encoded
        using a repetition code.

        Args:
            d (int): Number of code qubits (and hence repetitions) used.
            T (int): Number of rounds of ancilla-assisted syndrome measurement.


        Additional information:
            No measurements are added to the circuit if `T=0`. Otherwise
            `T` rounds are added, followed by measurement of the code
            qubits (corresponding to a logical measurement and final
            syndrome measurement round).
        """

        self.d = d
        self.mode = mode
        self.T = 0

        self.data_qubit = QuantumRegister(self.d**2,'data_qubit')
        self.ancilla_qubit = QuantumRegister((self.d**2-1),'ancilla_qubit')
        self.qubit_registers = {'data_qubit','ancilla_qubit'}

        self.ancilla_bits = []
        self.code_patches=[]
        self.side_code_patches=[]
        
        self.x_block_patches=[]
        self.z_block_patches=[]
        self.x_side_patches=[]
        self.z_side_patches=[]
        
        
        self.data_bit = ClassicalRegister(self.d**2,'code_bit')

        self.circuit = {}
        for log in ['0','1']:
            self.circuit[log] = QuantumCircuit(self.ancilla_qubit,self.data_qubit,name=log)

        self._preparation()
        
        for _ in range(T-1):
            self.stabilizer_circuit(reset=True)
            self.measurement()

        if T!=0:
            self.stabilizer_circuit(reset=False)
            self.measurement()
            
    def data_x(self, logs=('0', '1'), barrier=False):
        """
        Applies a logical x to the circuits for the given logical values.

        Args:
            logs (list or tuple): List or tuple of logical values expressed as
                strings.
            barrier (bool): Boolean denoting whether to include a barrier at
                the end.
        """
        for log in logs:
            for j in range(self.d):
                self.circuit[log].x(self.data_qubit[j])
            if barrier:
                self.circuit[log].barrier()
                
    def ancilla_H(self, logs=('0', '1'), barrier=False):
        """
        Applies a logical x to the circuits for the given logical values.

        Args:
            logs (list or tuple): List or tuple of logical values expressed as
                strings.
            barrier (bool): Boolean denoting whether to include a barrier at
                the end.
        """
        for log in logs:
            for j in range(self.d**2-1):
                self.circuit[log].h(self.ancilla_qubit[j])
            if barrier:
                self.circuit[log].barrier()
    
    def _preparation(self):
        """
        Prepares logical bit states by applying an x to the circuit that will
        encode a 1.
        """
        self.data_x(['1'])
        
    def code_patches_generation(self):
        """
        Generating the code patches for stabilizer circuit construction;
        
        Syndrome qubit assignment:
            0 to (d-1)^2-1 elements are the ones for the main block code patches;
            (d-1)^2 to d^2-1 elements are the ones for the side code patches, counted clockwise from the one above the first block code batch;
            
    
        """
        code_patches=[]
        code_patch_index = 0
        line_index = 0
        ancilla_index = 0
        
        # set up the ancilla_bits;
        self.ancilla_bits.append(ClassicalRegister((self.d**2-1),'round_'+str(self.T)+'_link_bit'))
        
        for j in range(self.d-1):
            for i in range(self.d-1):
                code_patch = {}
                qubit_index = line_index + i
                code_patch[code_patch_index] = [
                                                {0:self.data_qubit[qubit_index]}, {1:self.data_qubit[qubit_index+1]}, \
                                                {2:self.data_qubit[qubit_index+self.d]}, {3:self.data_qubit[qubit_index+self.d+1]}, \
                                                {4:self.ancilla_qubit[ancilla_index]}
                                               ]
                code_patches.append(code_patch)
                code_patch_index += 1
                ancilla_index += 1
            line_index += self.d
        
        x_block_patches=[]
        z_block_patches=[]
        
        if not (self.d%2):
            for i in range(len(code_patches)):
                if i%2:
                    x_block_patches.append(code_patches[i])
                else:
                    z_block_patches.append(code_patches[i])
        else:
            squared=[]
            for m in range(len(code_patches)//(self.d-1)):
                a=[]
                for k in range(self.d-1):
                    a.append(code_patches[k+(self.d-1)*m])
                squared.append(a)
        
            for m in range(self.d-1):
                if m%2:
                    for k in range(self.d-1):
                        if not (k%2):
                            x_block_patches.append(squared[m][k])
                        else:
                            z_block_patches.append(squared[m][k])
                else:
                    for k in range(self.d-1):
                        if not (k%2):
                            z_block_patches.append(squared[m][k])
                        else:
                            x_block_patches.append(squared[m][k])
                
        
        # create the side code batch index list;
        qubit_index_arr = [i for i in range(self.d**2)]
        side_patch_index = qubit_index_arr[:self.d]
        side_qubit_index = []
        for i in range(self.d-2):
            index=i*self.d+self.d
            side_index=[]
            for j in range(self.d):
                side_index.append(index+j)
            side_qubit_index.append(side_index[0])
            side_qubit_index.append(side_index[len(side_index)-1])
            index += self.d
        for i in range(len(side_qubit_index)):
            if i%2:
                side_patch_index.append(side_qubit_index[i])
        a=qubit_index_arr[len(qubit_index_arr)-self.d:]
        a.sort(reverse=True)
        side_patch_index = side_patch_index+a
        b=[]
        for i in range(len(side_qubit_index)):
            if (i+1)%2:
                b.append(side_qubit_index[i])
        b.sort(reverse=True)
        side_patch_index=side_patch_index+b
        
        # create the side patches;
        side_code_patches=[]
        side_code_patch_index = 0
        index=0
        for i in range((len(side_patch_index)//2)):
            code_patch = {}
            code_patch[side_code_patch_index]=[{0:self.data_qubit[side_patch_index[index]]},\
                                          {1:self.data_qubit[side_patch_index[index+1]]},\
                                          {2:self.ancilla_qubit[ancilla_index]}]
            index +=2
            side_code_patch_index +=1
            ancilla_index +=1
            side_code_patches.append(code_patch)
        
        top_num = self.d//2
        spacing = (len(side_code_patches) - 2*top_num)//2
        x_side_patches = side_code_patches[:top_num]+\
                            side_code_patches[len(side_code_patches)-spacing-top_num:len(side_code_patches)-spacing]
        
        z_side_patches = side_code_patches[top_num:top_num+spacing]+\
                            side_code_patches[top_num*2+spacing:]

        
        self.code_patches=code_patches
        self.side_code_patches=side_code_patches
        self.x_block_patches=x_block_patches
        self.x_side_patches=x_side_patches
        self.z_block_patches=z_block_patches
        self.z_side_patches=z_side_patches
        
        
    def stabilizer_circuit(self, reset=True):
        """
        Generating the circuits for the stabilizers of the surface code with NO measurements done;
        
        """
        self.code_patches_generation()
        
        # building the circuit;
        
        # initialize all the ancilla qubits for the x blocks;
        for log in ['0','1']:
            for i in range(len(self.x_block_patches)):
                key = [*self.x_block_patches[i]][0]
                ancilla = self.x_block_patches[i][key][len(self.x_block_patches[i][key])-1][len(self.x_block_patches[i][key])-1]
                self.circuit[log].h(ancilla)
            for i in range(len(self.x_side_patches)):
                key = [*self.x_side_patches[i]][0]
                ancilla = self.x_side_patches[i][key][len(self.x_side_patches[i][key])-1][len(self.x_side_patches[i][key])-1]
                self.circuit[log].h(ancilla)
                
        
        # the first cnot gate;
        
        # for all the block_patches;
        for log in ['0','1']:
            
            self.circuit[log].barrier()
            # for the x_blocks;
            for i in range(len(self.x_block_patches)):
                key = [*self.x_block_patches[i]][0]
                qubit = self.x_block_patches[i][key][0][0]
                ancilla = self.x_block_patches[i][key][len(self.x_block_patches[i][key])-1][len(self.x_block_patches[i][key])-1]
                self.circuit[log].cx(ancilla,qubit)
        
            
            # for the bottom x_sides;
            for i in range(len(self.x_side_patches[(self.d-1)//2:])):
                code_patch = self.x_side_patches[(self.d-1)//2:][i]
                key = [*code_patch][0]
                qubit = code_patch[key][1][1]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cx(ancilla,qubit)
            
            # for the z_blocks;
            for i in range(len(self.z_block_patches)):
                key = [*self.z_block_patches[i]][0]
                qubit = self.z_block_patches[i][key][0][0]
                ancilla = self.z_block_patches[i][key][len(self.z_block_patches[i][key])-1][len(self.z_block_patches[i][key])-1]
                self.circuit[log].cz(ancilla,qubit)
            
            # for the right_z_sides;
            for i in range(len(self.z_side_patches[:(self.d-1)//2])):
                code_patch = self.z_side_patches[:(self.d-1)//2][i]
                key = [*code_patch][0]
                qubit = code_patch[key][0][0]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cz(ancilla,qubit)
                
            #####################################
            #####################################
            
            self.circuit[log].barrier()
            # for the x_blocks;
            for i in range(len(self.x_block_patches)):
                key = [*self.x_block_patches[i]][0]
                qubit = self.x_block_patches[i][key][2][2]
                ancilla = self.x_block_patches[i][key][len(self.x_block_patches[i][key])-1][len(self.x_block_patches[i][key])-1]
                self.circuit[log].cx(ancilla,qubit)
        
            
            # for the top x_sides;
            for i in range(len(self.x_side_patches[:(self.d-1)//2])):
                code_patch = self.x_side_patches[:(self.d-1)//2][i]
                key = [*code_patch][0]
                qubit = code_patch[key][0][0]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cx(ancilla,qubit)
            
            # for the z_blocks;
            for i in range(len(self.z_block_patches)):
                key = [*self.z_block_patches[i]][0]
                qubit = self.z_block_patches[i][key][1][1]
                ancilla = self.z_block_patches[i][key][len(self.z_block_patches[i][key])-1][len(self.z_block_patches[i][key])-1]
                self.circuit[log].cz(ancilla,qubit)
            
            # for the left_z_sides;
            for i in range(len(self.z_side_patches[(self.d-1)//2:])):
                code_patch = self.z_side_patches[(self.d-1)//2:][i]
                key = [*code_patch][0]
                qubit = code_patch[key][1][1]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cz(ancilla,qubit)
        
            #####################################
            #####################################
            
            self.circuit[log].barrier()
            # for the x_blocks;
            for i in range(len(self.x_block_patches)):
                key = [*self.x_block_patches[i]][0]
                qubit = self.x_block_patches[i][key][1][1]
                ancilla = self.x_block_patches[i][key][len(self.x_block_patches[i][key])-1][len(self.x_block_patches[i][key])-1]
                self.circuit[log].cx(ancilla,qubit)
        
            
            # for the bottom x_sides;
            for i in range(len(self.x_side_patches[(self.d-1)//2:])):
                code_patch = self.x_side_patches[(self.d-1)//2:][i]
                key = [*code_patch][0]
                qubit = code_patch[key][0][0]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cx(ancilla,qubit)
            
            # for the z_blocks;
            for i in range(len(self.z_block_patches)):
                key = [*self.z_block_patches[i]][0]
                qubit = self.z_block_patches[i][key][2][2]
                ancilla = self.z_block_patches[i][key][len(self.z_block_patches[i][key])-1][len(self.z_block_patches[i][key])-1]
                self.circuit[log].cz(ancilla,qubit)
            
            # for the right_z_sides;
            for i in range(len(self.z_side_patches[:(self.d-1)//2])):
                code_patch = self.z_side_patches[:(self.d-1)//2][i]
                key = [*code_patch][0]
                qubit = code_patch[key][1][1]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cz(ancilla,qubit)
            
            #####################################
            #####################################
            
            self.circuit[log].barrier()
            # for the x_blocks;
            for i in range(len(self.x_block_patches)):
                key = [*self.x_block_patches[i]][0]
                qubit = self.x_block_patches[i][key][3][3]
                ancilla = self.x_block_patches[i][key][len(self.x_block_patches[i][key])-1][len(self.x_block_patches[i][key])-1]
                self.circuit[log].cx(ancilla,qubit)
        
            
            # for the top x_sides;
            for i in range(len(self.x_side_patches[:(self.d-1)//2])):
                code_patch = self.x_side_patches[:(self.d-1)//2][i]
                key = [*code_patch][0]
                qubit = code_patch[key][1][1]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cx(ancilla,qubit)
            
            # for the z_blocks;
            for i in range(len(self.z_block_patches)):
                key = [*self.z_block_patches[i]][0]
                qubit = self.z_block_patches[i][key][3][3]
                ancilla = self.z_block_patches[i][key][len(self.z_block_patches[i][key])-1][len(self.z_block_patches[i][key])-1]
                self.circuit[log].cz(ancilla,qubit)
            
            # for the left_z_sides;
            for i in range(len(self.z_side_patches[(self.d-1)//2:])):
                code_patch = self.z_side_patches[(self.d-1)//2:][i]
                key = [*code_patch][0]
                qubit = code_patch[key][0][0]
                ancilla = code_patch[key][len(code_patch[key])-1][len(code_patch[key])-1]
                self.circuit[log].cz(ancilla,qubit)
                
            for i in range(self.d**2-1):
                if reset:
                    self.circuit[log].reset(self.ancilla_qubit[i])
            self.circuit[log].barrier()
                
        self.T +=1
        
    def get_circuit_list(self):
        """
        Returns:
            circuit_list: self.circuit as a list, with
            circuit_list[0] = circuit['0']
            circuit_list[1] = circuit['1']
        """
        circuit_list = [self.circuit[log] for log in ['0', '1']]
        return circuit_list
    
            
    def measurement(self):
        """
        Measure the ancilla qubits of the circuit;
        """
        self.ancilla_bits.append(ClassicalRegister((self.d**2-1),'round_'+str(self.T)+'_link_bit'))
        for log in ['0','1']:
            self.circuit[log].add_register(self.ancilla_bits[-1])
            for i in range(self.d**2-1):
                self.circuit[log].measure(self.ancilla_qubit[i],self.ancilla_bits[self.T][i])