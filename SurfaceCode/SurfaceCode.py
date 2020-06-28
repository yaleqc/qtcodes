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
        self.data_bit = ClassicalRegister(self.d**2,'code_bit')

        self.circuit = {}
        for log in ['0','1']:
            self.circuit[log] = QuantumCircuit(self.ancilla_qubit,self.data_qubit,name=log)

        self._preparation()
        
        for _ in range(T-1):
            self.syndrome_measurement(reset=True,barrier=True)
        
        if T!=0:
            self.syndrome_measurement(reset=False)
            
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
            for j in range(self.d**2):
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
        self.ancilla_H()
        
    def syndrome_measurement(self, reset=True, barrier=False):
        # create the block code batches;
        code_patches=[]
        code_patch_index = 0
        line_index = 0
        ancilla_index = 0
        
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
            
        x_sequence = [0,1,2,3]
        z_sequence = [0,2,1,3]
        
        for log in ['0','1']:
            self.circuit[log].add_register(self.ancilla_bits[-1])
            for j in range(len(x_sequence)):
                for i in range(len(code_patches)):
                    ancilla = code_patches[i][i][len(code_patches[i][i])-1][len(code_patches[i][i])-1]
                    if i%2:
                        qubit = code_patches[i][i][x_sequence[j]][x_sequence[j]]
                        self.circuit[log].cx(ancilla,qubit)
                    else:
                        qubit = code_patches[i][i][z_sequence[j]][z_sequence[j]]
                        self.circuit[log].cz(ancilla,qubit)
    
            for i in range(len(side_code_patches)):
                index_1 = side_code_patches[i][i][0][0].index
                index_2 = side_code_patches[i][i][1][1].index
                if (index_1<self.d or index_1>(self.d**2-self.d)) and (index_2<self.d or index_2>(self.d**2-self.d-1)):
                    qubit_1 = side_code_patches[i][i][0][0]
                    qubit_2 = side_code_patches[i][i][1][1]
                    ancilla = side_code_patches[i][i][2][2]
                    self.circuit[log].cz(ancilla,qubit_1)
                    self.circuit[log].cz(ancilla,qubit_2)
                else:
                    qubit_1 = side_code_patches[i][i][0][0]
                    qubit_2 = side_code_patches[i][i][1][1]
                    ancilla = side_code_patches[i][i][2][2]
                    self.circuit[log].cx(ancilla,qubit_1)
                    self.circuit[log].cx(ancilla,qubit_2)
            self.circuit[log].barrier()
            for i in range(self.d**2-1):
                self.circuit[log].measure(self.ancilla_qubit[i],self.ancilla_bits[self.T][i])
                if reset:
                    self.circuit[log].reset(self.ancilla_qubit[i])
            self.circuit[log].barrier()
                
        self.T +=1
            
        
    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """
        for log in ['0', '1']:
            self.circuit[log].add_register(self.data_bit)
            self.circuit[log].measure(self.data_qubit, self.data_bit)