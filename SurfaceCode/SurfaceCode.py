# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 03:04:39 2020

@author: Andy
"""

from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit

class ModeInputException(Exception):
    def __init__(self, msg):
        super(ModeInputException, self).__init__(msg)

class SurfaceCode():
    """
    Implementation of a distance d repetition code, implemented over
    T syndrome measurement rounds.
    """

    def __init__(self, d, mode='normal'):
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

        self.data_qubit = QuantumRegister(self.d**2,'data_qubit')
        if mode=='normal':
            self.ancilla_qubit = QuantumRegister((self.d-1)**2,'ancilla_qubit')
        elif mode=='rotated':
            self.ancilla_qubit = QuantumRegister((self.d**2-1),'ancilla_qubit')
        else:
            raise ModeInputException('Input mode invalid. Should be either "normal" or "rotated".')
        self.qubit_registers = {'data_qubit','ancilla_qubit'}

        self.ancilla_bits = []
        self.data_bit = ClassicalRegister(self.d**2,'code_bit')

        self.circuit = {}
        for log in ['0','1']:
            self.circuit[log] = QuantumCircuit(self.ancilla_qubit,self.data_qubit,name=log)

        self._preparation()
        self.readout()
            
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
            if self.mode == 'normal':
                for j in range((self.d-1)**2):
                    self.circuit[log].h(self.ancilla_qubit[j])
                if barrier:
                    self.circuit[log].barrier()
            elif self.mode == 'rotated':
                for j in range(self.d**2-1):
                    self.circuit[log].h(self.ancilla_qubit[j])
                if barrier:
                    self.circuit[log].barrier()
            else:
                raise ModeInputException('Input mode invalid. Should be either "normal" or "rotated".')
    
    def _preparation(self):
        """
        Prepares logical bit states by applying an x to the circuit that will
        encode a 1.
        """
        self.data_x(['1'])
        self.ancilla_H()
        
    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """
        for log in ['0', '1']:
            self.circuit[log].add_register(self.data_bit)
            self.circuit[log].measure(self.data_qubit, self.data_bit)