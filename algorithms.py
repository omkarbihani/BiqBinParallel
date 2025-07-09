# TODO: docstings
# TODO: Add random seed preset for samplers to allow reproducibility.

# Required: Python 3.12

import logging
# import dwave.system
import numpy as np
import networkx as nx
import dwave
import dimod
import neal
import inspect
from abc import ABC, abstractmethod
import sys
import os
from copy import copy
from typing import Union, Tuple, Any
from collections import UserDict



logger = logging.getLogger()


type Graph = nx.Graph
type Sampler_response = Any 
type SimulatedAnnealingSampler = neal.SimulatedAnnealingSampler
type EmbeddingComposite = dwave.system.EmbeddingComposite
type Samplers = Union[SimulatedAnnealingSampler, EmbeddingComposite]
type AdjacencyMatrix = np.ndarray[tuple[int, int], np.dtype[np.int64]]
type QUBOMatrix = np.ndarray[tuple[int, int], np.dtype[np.float64]]


def get_A(filename: str) -> AdjacencyMatrix:
    """Read benchmark graph instance.

    Args:
        filename (str): file name.

    Returns:
        AdjacencyMatrix: Adjacency matrix of the graph.
    """
    A_raw = np.loadtxt(filename, dtype=int, usecols=(0, 1))
    n, k = A_raw[0, :]
    
    A = np.zeros(shape=(n,n))
    A[A_raw[1:,0]-1, A_raw[1:,1]-1] = 1
    A = A+A.transpose()

    return A_raw, A, k


class PrettyDict(UserDict):

  def __setitem__(self, key, value):
    if isinstance(value, dict):
        super().__setitem__(key, PrettyDict(value))
    else:        
        super().__setitem__(key, str(value))


class LazyEvaluationAsString:
    """Class for lazy evaluation of function as string
    """
    def __init__(self, func: callable, *args, **kwargs):
        """Initialize the class.

        Args:
            func (callable): A function.
        """
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        """Return string representation of func output.
        """
        return(f'{self.func(*self.args, **self.kwargs)}')


class PrintInfo:
    """
    Base class for info/debug printing.
    """
    class _LazyStrings:
        """Helper class for lazy joining multiple strings. 
        """
        def __init__(self, *args):
            self.args = args

        def __str__(self):
            return(''.join((str(s) for s in self.args)))
        
    def who(self, back:int = 0 ) -> str:
        """Return caller information.

        Args:
            back (int, optional): Caller frame id. Defaults to 0.

        Returns:
            str: Caller information. 
        """ 
        frame = sys._getframe( back + 1 )
        
        return f'{os.path.basename( frame.f_code.co_filename )}, {frame.f_lineno}, {type(self).__name__}.{frame.f_code.co_name}()'

    def print_info(self, message:str = '', *args):
        """Print info message.

        Args:
            message (str): Prints info message.
        """
        logging.info(
            self._LazyStrings(
                f'{self.who(1)}: ', # TODO: make it lazy ...
                message,
                *args
            )
        )
    
    def print_debug(self, message:str = '', *args):
        """Print debug message.

        Args:
            message (str): Prints debug message.
        """ 
        logging.debug(
            self._LazyStrings(
                f'{self.who(1)}: ', # TODO: make it lazy ...
                message,
                *args
            )
        )

    def print_warning(self, message:str = '', *args):
        """Print warning message.

        Args:
            message (str): Prints warning message.
        """
        logging.warning(
            self._LazyStrings(
                f'{self.who(1)}: ', # TODO: make it lazy ...
                message,
                *args
            )
        )
    
    def print_error(self, message:str = '', *args):
        """Print error message.

        Args:
            message (str): Prints debug error.
        """ 
        logging.error(
            self._LazyStrings(
                f'{self.who(1)}: ', # TODO: make it lazy ...
                message,
                *args
            )
        )


class Result:
    """Class to store results."""
    def __init__(self, sampler_response: Sampler_response, g: Graph, algorithm: str, **additional_info):
        self._data = {'g': g.copy(), 'sampler_response': sampler_response, 'additional_info': copy(additional_info), 'algorithm': str(algorithm)}

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        return self._data
    
    def __getitem__(self, item):
         return self._data[item]
    
    def __setitem__(self, key, value): 
        self._data[key] = value
    
    def update(self, data:dict):
        self._data = {**self._data, **data}

    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return (
            f'G={self._data['g'].number_of_nodes(), self._data['g'].number_of_edges()}\n'
            f'f(x)={self._data['sampler_response']['value']}\n'
            f'algorithm={self._data['algorithm']}'
        )
    

class QuboTriU(PrintInfo):
    """Class for QUBO representation in an upper triangular form.
    """
    def __init__(self, g: Graph):
        """Initialize the class.

        Args:
            g (Graph): Graph.
        """
        self.print_debug()
        
        self.A = nx.adjacency_matrix(g).todense().astype(float)

    def update(self, l: float, mu: float=0, k: int=0) -> QUBOMatrix:
        """Update QUBO matrix

        Args:
            l (float): Lambda parameter.
            mu (float, optional): Mu parameter. Defaults to 0.
            k (int, optional): Target cardinality of the sparsest subgraph.. Defaults to 0.

        Returns:
            QUBOMatrix: Updated QUBO matrix.
        """
        self.print_debug()

        # !!!! triangle structure is better performing!!!
        Q = self.A + mu
        Q = np.triu(Q, 1)
        np.fill_diagonal(Q, -l+mu/2-mu*k)
        self.Q = Q

        return self.Q 

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        return self
    
    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return f'(QuboTriU: Qubo in upper triangular form, A={self.A.shape})'
    

class QuboSym(QuboTriU):
    """Class for representing QUBO as a symmetric matrix.
    """
    def __init__(self, g: Graph): 
        """Initialize the class.

        Args:
            g (Graph): Graph.
        """
        super().__init__(g)
        
        self.A = nx.adjacency_matrix(g).todense().astype(float)
    
    def update(self, l: float, mu: float=0, k: int=0) -> QUBOMatrix:
        """Update QUBO matrix

        Args:
            l (float): Lambda parameter.
            mu (float, optional): Mu parameter. Defaults to 0.
            k (int, optional): Target cardinality of the sparsest subgraph.. Defaults to 0.

        Returns:
            QUBOMatrix: Updated QUBO matrix.
        """
        super().update(l, mu, k)
        self.Q = 1/2*(self.Q+self.Q.T)

        return self.Q

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        return self
    
    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return f'(QuboSym: Qubo as symmetric matrix, A={self.A.shape})'


type QUBO = Union[QuboTriU, QuboSym]


class Sampler(PrintInfo):
    """Compute approximation of a minimization problem defined by QUBO.
    """
    def __init__(self, sampler: Samplers, num_reads: int):
        self.print_debug()

        self.sampler = sampler
        self.num_reads = num_reads
    def compute(self, Q: QUBOMatrix) -> Sampler_response:
        """_summary_

        Args:
            Q (QUBOMatrix): QUBO matrix,

        Returns:
            Sampler_response: Sampler response + additional information.
        """
        self.print_debug()

        bqm = dimod.BQM.from_qubo(Q)
        sampleset = self.sampler.sample(bqm, num_reads = self.num_reads)
        
        X = sampleset.first.sample
        value = sampleset.first.energy
        x = np.zeros(len(X), dtype=int)
        for i, v in X.items():
            x[i] = v
        return {'x': x, 'value': value, 'cardinality': sum(x), 'sampler_response': sampleset}

    def set_num_reads(self, num_reads: int):
        """Set num_reads parameter for the sampler.

        Args:
            num_reads (int): Number of reads.
        """
        self.num_reads = num_reads

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        return {
            "sampler": str(self)
        }
    
    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return f'{type(self).__name__}(sampler={type(self.sampler).__name__}(), num_reads={self.num_reads})'


class DidNotConverge(Exception):
    """Exception raised if algorithm did not converge."""

    def __init__(self, message):
        self.message = message
        super().__init__(message)

    def __str__(self):
        return f'Algorithm did not converge. Try to adjust input parameters.\n{self.message}'    
    

class Algorithm(ABC, PrintInfo):
    """Abstract base class for k-sparsest subgraph algorithm. 
    """
    def __init__(self, g: Graph, method: Sampler, result: Result, qubo: QUBO=QuboTriU):
        """Initialize the class.

        Args:
            g (Graph): Input graph instance.
            method (Sampler): Sampler to use for computation.
            result (Result): Class to be used to sore result.
            qubo (QUBO, optional): Class for QUBO representation. Defaults to QuboTriU.
        """
        self.print_debug()

        # Convert node labeling to integers starting from 0. Store mapping to new labeling in '_original_labeling' node attribute.
        # Save a copy of the graph.
        self.g = nx.convert_node_labels_to_integers(g, first_label=0, ordering='default', label_attribute='_original_labeling')
        self.method = method
        self.qubo = qubo(self.g)
        self.result = result
        self.computed_cardinality = None
        self.mui = 0
        self.best_computed_cardinality = None
        self.best_computed_cardinality_x = None
        self.best_sampler_response = None
        self.k = None
   
    def compute(self) -> Result:
        """Compute an approximation for the k-sparsest subgraph problem.

        Returns:
            Result: Result.
        """
        while(self.computed_cardinality!=self.k and self.update_lambda_mu(self.computed_cardinality)):
            Q = self.qubo.update(l=self.lambdai, mu=self.mui, k=self.k)
            sampler_response = self.method.compute(Q)
            self.computed_cardinality = sampler_response['cardinality']  
            
            self.print_debug(f'k = {self.k}, computed_k = {self.computed_cardinality}, lam = {self.lambdai}, mu = {self.mui}')

            if (self.best_computed_cardinality is None) and (self.computed_cardinality >= self.k):
                self.best_computed_cardinality = self.computed_cardinality
                self.best_computed_cardinality_x = sampler_response['x']
                self.best_sampler_response = sampler_response

            if (self.computed_cardinality >= self.k) and (self.computed_cardinality < self.best_computed_cardinality):
                self.best_computed_cardinality = self.computed_cardinality
                self.best_computed_cardinality_x = sampler_response['x']
                self.best_sampler_response = sampler_response

        if self.best_computed_cardinality is None:
            raise DidNotConverge(self.did_not_converge_message())

        g = nx.induced_subgraph(self.g, np.argwhere(self.best_computed_cardinality_x==1).flatten()).copy()
        
        return self.result(sampler_response=self.best_sampler_response, g=g, algorithm=str(self))

    @abstractmethod
    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        return {
            'algorithm': f'{type(self).__name__}():',
            'g': self.g,
            'method': self.method.summary(),
            'qubo': self.qubo.summary(),
            'computed_cardinality': self.computed_cardinality,
            'mui': self.mui,
            'best_computed_cardinality': self.best_computed_cardinality,
            'best_computed_cardinality_x': self.best_computed_cardinality_x,
        }
       
    @abstractmethod
    def update_lambda_mu(self, k: int) -> bool:
        """Update lambda and mu parameters.

        Args:
            k (int): Size of k-sparsest subgraph.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        ...

    def did_not_converge_message(self) -> str:
        return str(self) 

    @abstractmethod
    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
        f'{type(self).__name__}'
        f'(g={self.g.number_of_nodes(), self.g.number_of_edges()},'
        f'k={self.k}, method={self.method}, result={self.result.__name__}, qubo={self.qubo}'
        )

class AugumentedLagrangian(Algorithm):
    """Class implementing the Augemented lagrangian method for the k-sparsest subgraph problem.
    """
    def __init__(self, g: Graph, k: int, method: Sampler, result: Result, 
                 initial_lambda: float, initial_mu: float, step: float, max_iterations: int):
        """Initialize the class.

        Args:
            g (Graph): Input graph instance.
            k (int): Target cardinality of the sparsest subgraph.
            method (Sampler): Sampler method to use for computation.
            result (Result): Class for storing result.
            initial_lambda (float): Initial lambda parameter value.
            initial_mu (float): Initial mu parameter value. Can't be 0.
            step (float): Step value.
            max_iterations (int): Maximumi terations.
        """
        super().__init__(g, method, result)
        self.k = k
        self.lambdai = initial_lambda
        self.initial_lambda = initial_lambda
        
        if initial_mu <= 0:
            raise ValueError("mu cannot be less than or equal to 0.") 
        
        self.mui = initial_mu

        self.initial_mu = initial_mu
        self.max_iterations = max_iterations
        self.number_of_iterations = 0
        self.step = step

    
    def update_lambda_mu(self, k: int) -> bool:
        """Update lambda and mu parameters.

        Args:
            k (int): Size of k-sparsest subgraph.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        self.print_debug(f"iter = {self.number_of_iterations}")
        if k is None:
            self.number_of_iterations += 1

            return True
        
        if self.number_of_iterations < self.max_iterations:
            self.lambdai += self.mui*(self.k-k)
            self.mui *= self.step
            self.number_of_iterations += 1
            
            return True
        return False           

    def compute(self) -> Result:
        """Compute an approximation for the k-sparsest subgraph problem using the Augemented langrangian method.

        Returns:
            Result: Result.
        """
        result = super().compute()
        
        if self.best_computed_cardinality != self.k:
            raise DidNotConverge(self.did_not_converge_message())
        
        result['additional_info']['algo_summary'] = self.summary()
        
        return result    

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: _description_summary
        """
        summary = super().summary()
        summary['initial_lambda'] = self.initial_lambda
        summary['lambdai'] = self.lambdai
        summary['mui'] = self.mui
        summary['step'] = self.step
        summary['max_iterations'] = self.max_iterations
        summary['number_of_iterations'] = self.number_of_iterations

        return summary

    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
            f'{type(self).__name__}'
            f'(g={self.g.number_of_nodes(), self.g.number_of_edges()},'
            f'k={self.k}, method={self.method}, result={self.result.__name__}, initial_lambda={self.initial_lambda},'
            f'initial_mu={self.initial_mu}, step={self.step}, max_iterations={self.max_iterations}))\n'
            f'lambdai={self.lambdai}\n'
            f'mui={self.mui}\n'
            f'number_of_iterations={self.number_of_iterations}'
        )


class Linear(Algorithm):
    """Class implementing the Linear version of algorithm for k-sparsest subgraph problem.
    If target cardinality k is not reached, the greedy algorithm is used from the closest upper k found.

    Args:
        Algorithm (Algorithm): Inherits from the base class Algorithm.
    """
    def __init__(self, g: Graph, k: int, method: Sampler, result: Result, 
                 initial_lambda: float, lambda_max_iterations: int, lambda_step: float):
        """Iitialize the class.

        Args:
            g (Graph): Graph instance.
            k (int): Target cardinality of k-sparsest subgraph.
            method (Sampler): Sampler to use.
            result (Result): Class to store result.
            initial_lambda (float): Initital lambda value.
            lambda_max_iterations (int): Maximum number of iterations.
            lambda_step (float): Lambda step value.
        """
        super().__init__(g, method, result)
        self.print_debug()

        self.k = k
        self.lambdai = initial_lambda
        self.initial_lambda = initial_lambda
        self.lambda_max_iterations = lambda_max_iterations
        self.lambda_number_of_iterations = 0
        self.lambda_step = lambda_step
        self.best_computed_cardinality = None
        self.best_computed_cardinality_x = None
        self.greedy_used = False
    
    def update_lambda_mu(self, k: int) -> bool:
        """Update lambda and mu parameters.

        Args:
            k (int): Size of k-sparsest subgraph.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        self.print_debug()

        if k is None:
            self.lambda_number_of_iterations += 1
            return True
        if self.lambda_number_of_iterations < self.lambda_max_iterations:
            self.lambdai += self.lambda_step*(self.k-k)
            self.lambda_number_of_iterations += 1
            return True
        return False

    def compute(self) -> Result:
        """Find an approximation of the sparsest subgraph using lagrangian.

        Returns:
            Result: Result.
        """
        result = super().compute()
        g = result['g']
        
        if self.best_computed_cardinality != self.k:
            self.greedy_used = True
            self.print_info(f'Didn\'t find target cardinality={self.k}, Going greedy with k={self.best_computed_cardinality}')
            G = nx.induced_subgraph(self.g, np.argwhere(self.best_computed_cardinality_x==1).flatten())
            g_greedy, dmax, davg = greedy(G, self.k)
            result['g'] = g_greedy
            result['additional_info']['g_for_greedy'] = g
            # Update pretty algorithm representaion. We changed some state parameters.
            # TODO: do it be better ....
            result['algorithm'] = str(self)
        
        result['additional_info']['algo_summary'] = self.summary()
               
        return result

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        summary = super().summary()
        summary['initial_lambda'] = self.initial_lambda
        summary['lambdai'] = self.lambdai
        summary['lambda_max_iterations'] = self.lambda_max_iterations
        summary['lambda_number_of_iterations'] = self.lambda_number_of_iterations
        summary['lambda_step'] = self.lambda_step
        summary['best_computed_cardinality'] = self.best_computed_cardinality
        summary['greedy_used'] = self.greedy_used

        return summary

    def __str__(self) -> str:
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
            f'{type(self).__name__}'
            f'(g={self.g.number_of_nodes(), self.g.number_of_edges()},'
            f'k={self.k}, method={self.method}, result={self.result.__name__}, initial_lambda={self.initial_lambda},'
            f'lambda_max_iterations={self.lambda_max_iterations}, lambda_step={self.lambda_step}))\n'
            f'lambdai={self.lambdai}\n'
            f'greedy_used={self.greedy_used}\n'
            f'number_of_iterations={self.lambda_number_of_iterations}\n'
            f'best_computed_cardinality={self.best_computed_cardinality}'
        )


class ProposedAlgorithm(Algorithm):
    """ 
    Implements a proposed approach to find sparsest-k-subgraph 
    
    This method first searches for an optimal lambda using bisection, and if it fails 
    to find the required subgraph, it increases the quadratic penalty term (mu) iteratively.
    """

    def __init__(self, g: Graph, k: int, method: Sampler, result: Result, 
                 initial_lambda: float, initial_lambda_min: float, initial_lambda_max:float, lambda_max_iterations: int,
                 mu_step: float, mu_max_iterations: int, mu_num_reads: int):
        super().__init__(g, method, result)
        self.print_debug()

        self.k = k
        self.initial_lambda = initial_lambda

        if initial_lambda_min > initial_lambda_max:
            raise ValueError("lambda_min must be greater than lambda_max")

        self.initial_lambda_min = initial_lambda_min
        self.initial_lambda_max = initial_lambda_max
        self.lambda_min = initial_lambda_min
        self.lambda_max = initial_lambda_max
        self.lambda_max_iterations = lambda_max_iterations
        self.mu_step = mu_step
        self.mu_max_iterations = mu_max_iterations
        self.mu_num_reads = mu_num_reads

        self.lambda_number_of_iterations = 0
        self.mu_number_of_iterations = 0
        self.lambdai = initial_lambda

    def update_lambda_mu(self, k):
        """
         Updates the values of lambda and mu based on the current subgraph size.
        - If bisection search on lambda fails to find the correct k, mu is increased iteratively.

        Args:
            k (int): Size of k-sparsest subgraph.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        self.print_debug(f"lambda_iter = {self.lambda_number_of_iterations}, mu_iter={self.mu_number_of_iterations}, lambdai = {self.lambdai}, mui=  {self.mui}")

        if (self.lambda_number_of_iterations < self.lambda_max_iterations):
            
            if k is None:
                self.lambda_number_of_iterations += 1
                return True
            if k > self.k:
                self.lambda_max = self.lambdai
            if k < self.k:
                self.lambda_min = self.lambdai
            if k != self.k:
                self.lambdai = (self.lambda_min + self.lambda_max)/2
            
            self.lambda_number_of_iterations += 1

            return True
        
        elif self.mu_number_of_iterations < self.mu_max_iterations:
            self.method.set_num_reads(self.mu_num_reads)
            self.mui += self.mu_step
            self.mu_number_of_iterations += 1

            return True
        else:
            return False
        
    def compute(self) -> Result:
        """Compute an approximation for the k-sparsest subgraph problem using the Proposed method.

        Returns:
            Result: Result.
        """
        result = super().compute()
        
        if self.best_computed_cardinality != self.k:
            raise DidNotConverge(self.did_not_converge_message())
            
        result['additional_info']['algo_summary'] = self.summary()
        return result                
        
    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: _description_summary
        """
        summary = super().summary()
        summary['initial_lambda'] = self.initial_lambda
        summary['lambdai'] = self.lambdai
        summary['lambda_number_of_iterations'] = self.lambda_number_of_iterations
        summary['lambda_max_iterations'] = self.lambda_max_iterations
        summary['mui'] = self.mui
        summary['mu_step'] = self.mu_step
        summary['mu_number_of_iterations'] = self.mu_number_of_iterations 
        summary['mu_max_iterations'] = self.mu_max_iterations
        summary['mu_num_reads'] = self.mu_num_reads
        summary['initial_lambda_max'] = self.initial_lambda_max
        summary['initial_lambda_min'] = self.initial_lambda_min
        summary['lambda_max'] = self.lambda_max
        summary['lambda_min'] = self.lambda_min

        return summary

    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
            f'{type(self).__name__}'
            f'(g={self.g.number_of_nodes(), self.g.number_of_edges()}), '
            f'k={self.k}, method={self.method}, result={self.result.__name__}, \n'
            f'initial_lambda={self.initial_lambda}, initial_lambda_min={self.initial_lambda_min}, '
            f'initial_lambda_max={self.initial_lambda_max}, lambda_max_iterations={self.lambda_max_iterations} '
            f'mu_step={self.mu_step}, mu_max_iterations = {self.mu_max_iterations}, mu_num_reads={self.mu_num_reads}))\n' 
            f'lambdai={self.lambdai}, mui={self.mui}\n'
            f'lambda_number_of_iterations={self.lambda_number_of_iterations}, mu_number_of_iterations={self.mu_number_of_iterations}\n'
            f'lambda_min={self.lambda_min}, lambda_max={self.lambda_max}, best_computed_cardinality={self.best_computed_cardinality}'
        )

class QuadraticPenalty(Algorithm):
    def __init__(self, g: Graph, k: int, initial_mu: float, step: float, mu_max_iterations: int, method: Sampler):
        super().__init__(g=g, method=method, result=Result)
        self.initial_mu = initial_mu
        self.mui = initial_mu
        self.lambdai = 0
        self.step = step
        self.k = k
        self.mu_max_iterations = mu_max_iterations
        self.mu_number_of_iterations = 0
    

    def compute(self) -> Result:
        """Compute an approximation for the k-sparsest subgraph problem using the Augemented langrangian method.

        Returns:
            Result: Result.
        """
        result = super().compute()
        
        if self.best_computed_cardinality != self.k:
            raise DidNotConverge(self.did_not_converge_message())
        
        result['additional_info']['algo_summary'] = self.summary()
        
        return result 
    
    def update_lambda_mu(self, k: int) -> bool:
        """Update lambda and mu parameters.

        Args:
            k (int): Size of k-sparsest subgraph.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        self.print_debug(f"iter = {self.mu_number_of_iterations}")
        if k is None:
            self.mu_number_of_iterations += 1

            return True
        
        if self.mu_number_of_iterations < self.mu_max_iterations:
            self.mui += self.step
            self.mu_number_of_iterations += 1
            
            return True
        return False  

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: _description_summary
        """
        summary = super().summary()
        summary['initial_mu'] = self.initial_mu
        summary['mui'] = self.mui
        summary['mu_number_of_iterations'] = self.mu_number_of_iterations
        summary['mu_max_iterations'] = self.mu_max_iterations
        return summary

    def __str__(self):
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
            f'{type(self).__name__}'
            f'(g={self.g.number_of_nodes(), self.g.number_of_edges()},'
            f'k={self.k}, method={self.method}, result={self.result.__name__},'
            f'initial_mu={self.initial_mu}, step={self.step}, mu_max_iterations={self.mu_max_iterations}))\n'
            f'lambdai={self.lambdai}, mui={self.mui} \n'
            f'mu_number_of_iterations={self.mu_number_of_iterations}'
        )
               
class MaxIndependentSet(Algorithm):
    """Class implementing the Independent set of the graph.

    Args:
        Algorithm (Algorithm): Inherits from the base class Algorithm.
    """
    def __init__(self, g: Graph, method: Sampler, result: Result, initital_lambda: float=None,
                 lambda_min: float=0 , lambda_max: float=1,  lambda_max_iterations: int=10):
        """Iitialize the class.

        Args:
            g (Graph): Graph instance.
            method (Sampler): Sampler to use.
            result (Result): Class to store result.
            lambda_min (float): Minimum lambda value. Default = 0
            lambda_max (float): Maximum lambda value. Default = 1
            lambda_max_iterations (int): Maximum number of iterations. Default=10
        """
        super().__init__(g, method, result)
        self.print_debug()

        if (lambda_min < 0 or lambda_max > 1 or lambda_min > lambda_max):
            raise ValueError("For MIS problem, please set lambda_min>=0, lambda_max<=1 and lambda_min<=lambda_max")
        
        self.initial_lambda = initital_lambda
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lambda_max_iterations = lambda_max_iterations
        self.lambda_number_of_iterations = 0
        self.best_computed_cardinality = None
        self.best_computed_cardinality_x = None
        self.networkx_heuristic_used = False


    def calculate_edges(self, sampler_response):
        return 0.5*sampler_response['x'] @ self.qubo.A @ sampler_response['x']
    
    def update_lambda_mu(self) -> bool:
        """
        Update lambda and mu parameters.
        
        This method randomly adjusts lambda within the specified bounds (`lambda_min` and `lambda_max`) 
        to try and converge towards an optimal independent set. It stops adjusting lambda after the maximum 
        allowed iterations.

        Returns:
            bool: Returns True if updated, False otherwise (if max iteration reached, for example). 
        """
        if self.lambda_number_of_iterations < self.lambda_max_iterations:
            if self.initial_lambda is not None and self.lambda_number_of_iterations == 0:
                self.lambdai = self.initial_lambda
            else:
                self.lambdai = np.random.uniform(low=self.lambda_min, high=self.lambda_max)
            self.lambda_number_of_iterations += 1
            return True
        return False

    def compute(self) -> Result:
        while(self.update_lambda_mu()):
            Q = self.qubo.update(l=self.lambdai, mu=0, k=0)
            sampler_response = self.method.compute(Q)
            self.computed_cardinality = sampler_response['cardinality']  
            
            self.print_debug(f"computed_k = {self.computed_cardinality}, lam = {self.lambdai} num_edges = ", 
                             #  This a neat trick to lazy evaluate expensive debug computations.
                             #  It would be evaluated in non DEBUG modes otherwise.
                             LazyEvaluationAsString(self.calculate_edges, sampler_response))

            if (self.best_computed_cardinality is None) or (self.best_computed_cardinality < self.computed_cardinality):
                self.best_computed_cardinality = self.computed_cardinality
                self.best_computed_cardinality_x = sampler_response['x']
                self.best_sampler_response = sampler_response
        
        g = nx.induced_subgraph(self.g, np.argwhere(self.best_computed_cardinality_x==1).flatten()).copy()

        result = self.result(sampler_response=self.best_sampler_response, g=g, algorithm=str(self))

        if g.number_of_edges() != 0:
            self.print_info("Couldn't find the independent set. Computing it by using networkx algorithm on the best found induced subgraph")
            ids_set = nx.algorithms.approximation.maximum_independent_set(g)
            g_nx = nx.induced_subgraph(g,ids_set)
            self.networkx_heuristic_used = True

            result['additional_info']['g_for_nx'] = g
            result['g'] = g_nx
            # Update pretty algorithm representaion. We changed some state parameters.
            # TODO: do it be better ....
            result['algorithm'] = str(self)
        
        result['additional_info']['algo_summary'] = self.summary()

        return result
    

    def summary(self) -> dict:
        """Return summary of stored information.

        Returns:
            dict: Summary.
        """
        summary = super().summary()
        summary['lambdai'] = self.lambdai
        summary['lambda_min'] = self.lambda_min
        summary['lambda_max'] = self.lambda_max
        summary['lambda_max_iterations'] = self.lambda_max_iterations
        summary['lambda_number_of_iterations'] = self.lambda_number_of_iterations
        summary['initial_lambda'] = self.initial_lambda
        summary['networkx_heuristic_used'] = self.networkx_heuristic_used

        return summary

    def __str__(self) -> str:
        """Represent the class as pretty string.

        Returns:
            str: Pretty string representation of the class.
        """
        return ( 
            f'{type(self).__name__}'
            f'(g={self.g.number_of_nodes(), self.g.number_of_edges()},'
            f'best_computed_cardinality={self.best_computed_cardinality}, method={self.method}, '
            f'result={self.result.__name__}, initial_lambda={self.initial_lambda}, lambda_min={self.lambda_min},'
            f'lambda_max={self.lambda_max}, lambda_max_iterations={self.lambda_max_iterations}))\n'
            f'Qubo={self.qubo}\n'
            f'lambdai={self.lambdai}, lambda_number_of_iterations={self.lambda_number_of_iterations}, networkx_heuristic_used={self.networkx_heuristic_used}'
        )


def greedy(g: Graph, k: int) -> Tuple[Graph, int, float]:
    """Find k-sparsest subgraph using greedy algorithm.

    Args:
        g (Graph): Graph.
        k (int): Size of a sparsest subrgraph.

    Returns:
        Tuple[Graph, int, float]: (Approximation of k-sparsest subgraph, maximum degree, average degree)
    """
    G = g.copy()
    n = G.number_of_nodes()

    for i in range(n-k):
        degree_sequence = sorted(((d, n) for n, d in G.degree()), reverse=True)
        dmax, node = degree_sequence[0]
        G.remove_node(node)

    degree_sequence = sorted(((d, n) for n, d in G.degree()), reverse=True)
    dmax, node = degree_sequence[0]
    davg = G.number_of_edges()/G.number_of_nodes()/2
    
    return G, dmax, davg

def relabel_subgraph_nodes(g: Graph, labeling_map: dict) -> Graph:
    """Relabel nodes of g using the map: labeling_map

    Args:
        g (Graph): A graph.
        labeling_map (dict): Node labeling map.

    Returns:
        Graph: Relabeled graph g.
    """
    g_node_labels = g.nodes
    try:
        mapping = {label: labeling_map[label] for label in g_node_labels}
    except KeyError:
        logger.error('Some node labels of graph g are not in labeling_map.')
        logger.error(f'labelig_map={labeling_map}')
        logger.error(f'g_node_labels={g_node_labels}')
        raise 
    except Exception as e:
        logger.error(f'Unexpected {type(e)}')
        raise

    return nx.relabel_nodes(g, mapping, copy=True)
