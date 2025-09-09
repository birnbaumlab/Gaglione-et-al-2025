import numpy as np
import random
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import Parallel, delayed
from Bio.Seq import Seq
import os
import json



class OverhangSet:
    """
    Class to hold a set of overhangs for a TCRAFT cloning scheme.
    Can be mutable (for optimization) or immutable (for final output).
    Supports random initialization, mutation, and indexing by variable ID or index.
    Also supports saving/loading to/from dictionary or CSV/JSON files.
    4bp overhangs only, no palindromic or self-complementary overhangs allowed.
    Must specify whether V is on left or right of C in the cloning scheme.
    Must specify restriction enzyme (BsaI, BbsI, or BsmBI) for scoring purposes.
    """

    order_VLeft_CRight = "V on Left, C on Right"
    order_CLeft_VRight = "C on Left, V on Right"
    
    # Constructors for OverhangSet-------

    def __init__(self, v_seqs, c_seq, v_ids, c_id, restriction_enzyme, order, is_mutable=True, random_init=True, banned_OHs=[]):
        """
        Initialize an OverhangSet object.
        v_seqs: list of variable region protein sequences (strings) - should be last 7-8 amino acids of V region
        c_seq: constant region protein sequence (string) - should be first 8-10 amino acids of C region
        v_ids: list of variable region IDs (strings) - should match length of v_seqs
        c_id: constant region ID (string)
        restriction_enzyme: string, one of 'BsaI', 'BbsI', or 'BsmBI'
        order: string, either OverhangSet.order_VLeft_CRight or OverhangSet.order_CLeft_VRight
        is_mutable: bool, whether the OverhangSet is mutable (for optimization) or immutable (final output)
        random_init: bool, if True and is_mutable, will randomly initialize the overhangs
        banned_OHs: list of overhangs (strings) that should not be used
        """

        self.v_seqs = v_seqs
        self.c_seq = c_seq
        self.v_ids = v_ids
        self.c_id = c_id
        self.restriction_enzyme = restriction_enzyme
        if order not in [OverhangSet.order_VLeft_CRight, OverhangSet.order_CLeft_VRight]:
            raise ValueError("Invalid order specified. Must use 'OverhangSet.order_VLeft_CRight' or 'OverhangSet.order_CLeft_VRight'.")
        self.order = order
        self.is_mutable = is_mutable
        self.banned_OHs = banned_OHs
        
        if self.is_mutable:
            self._prepare_valid_overhangs()
            if random_init:
                self.random_initialize()

    
    @classmethod
    def build(cls, full_v_seqs, full_c_seq, v_ids, c_id, restriction_enzyme, order, v_window=7, c_window=9, banned_OHs=[]):
        """
        Build an OverhangSet from full protein sequences by extracting the last v_window amino acids of each V sequence
        and the first c_window amino acids of the C sequence.
        full_v_seqs: list of full variable region protein sequences (strings)
        full_c_seq: full constant region protein sequence (string)
        v_ids: list of variable region IDs (strings) - should match length of v_seqs
        c_id: constant region ID (string)
        restriction_enzyme: string, one of 'BsaI', 'BbsI', or 'BsmBI'
        order: string, either OverhangSet.order_VLeft_CRight or OverhangSet.order_CLeft_VRight
        v_window: int, number of amino acids to take from end of V sequences
        c_window: int, number of amino acids to take from start of C sequence
        banned_OHs: list of overhangs (strings) that should not be used

        Returns an OverhangSet object constructed from these sequences.
        """
        v_seqs = [i[-v_window:] for i in full_v_seqs]
        c_seq = full_c_seq[:c_window]
        return cls(v_seqs, c_seq, v_ids, c_id, restriction_enzyme, order, is_mutable=True, random_init=True, banned_OHs=banned_OHs)


    @classmethod
    def create_from_existing_overhangs(cls, v_oh_list, c_oh_list, v_ids, c_id, restriction_enzyme, order, is_mutable=False, banned_OHs=[]):
        """
        Create an OverhangSet from existing overhang sets.
        v_oh_list: list of variable region 4bp overhangs (strings)
        c_oh_list: list of constant region 4bp overhangs (strings)
        v_ids: list of variable region IDs (strings) - should match length of v_seqs
        c_id: constant region ID (string)
        restriction_enzyme: string, one of 'BsaI', 'BbsI', or 'BsmBI'
        order: string, either OverhangSet.order_VLeft_CRight or OverhangSet.order_CLeft_VRight
        is_mutable: bool, whether the OverhangSet is mutable (for optimization) or immutable (final output)
        banned_OHs: list of overhangs (strings) that should not be used
        Returns an OverhangSet object constructed from these overhangs.
        """
        obj = cls([], '', v_ids, c_id, restriction_enzyme, order, is_mutable, random_init=False, banned_OHs=banned_OHs)
        obj.v_oh_set = np.array(v_oh_list, dtype='<U4')
        obj.c_oh_set = np.array(c_oh_list, dtype='<U4')
        return obj
    

    @classmethod
    def load_from_dict(cls, data):
        """
        Load an OverhangSet from a dictionary (as produced by dump_to_dict).
        data: dictionary with keys matching those produced by dump_to_dict
        Returns an OverhangSet object constructed from this dictionary.
        """

        if not isinstance(data, dict):
            raise ValueError("Input data must be a dictionary.")
        
        v_seqs = data.get('v_seqs', [])
        c_seq = data.get('c_seq', '')
        v_ids = data.get('v_ids', [])
        c_id = data.get('c_id', '')
        restriction_enzyme = data.get('restriction_enzyme', 'BsaI')
        order = data.get('order', cls.order_VLeft_CRight)
        is_mutable = data.get('is_mutable', True)
        banned_OHs = data.get('banned_OHs', [])
        
        obj = cls(v_seqs, c_seq, v_ids, c_id, restriction_enzyme, order, is_mutable, random_init=False, banned_OHs=banned_OHs)
        obj.v_oh_set = np.array(data.get('v_oh_set', []), dtype='<U4')
        obj.c_oh_set = np.array(data.get('c_oh_set', []), dtype='<U4')
        
        return obj
    
    # Main public methods for OverhangSet-------

    def random_initialize(self):
        """
        Randomly initialize the overhang set, ensuring no duplicates or self-complementary overhangs
        """

        if not self.is_mutable:
            raise ValueError("Cannot random initialize a non-mutable OverhangSet.")
        
        self.v_oh_set = np.empty((len(self.v_seqs),), dtype='<U4')
        self.c_oh_set = np.empty((len(self.v_seqs),), dtype='<U4')

        for i in range(self.v_oh_set.shape[0]):
            while True:
                #start by picking a constant region overhang
                c_choice = random.choice(self.valid_c_overhangs)
                if self._not_used_yet(c_choice):
                    self.c_oh_set[i] = c_choice
                    # once that worked, then pick a variable region overhang
                    completed = False
                    while True:
                        v_choice = random.choice(self.valid_v_overhangs[i])
                        if self._not_used_yet(v_choice):
                            self.v_oh_set[i] = v_choice
                            completed = True
                            break
                    if completed:
                        break
            

    def mutate(self, idx, side):
        """
        Mutate the overhang at position idx (int or V_ID string) on side ('V' or 'C').
        Returns True if mutation was successful, False if no valid mutation was possible.
        idx: int [0, len(v_ids)-1] or string (V_ID, C_ID)
        side: int [0, 1] or string ('V', 'C')
        Returns True if mutation was successful, False if no valid mutation was possible.
        """
        if not self.is_mutable:
            raise ValueError("Cannot mutate a non-mutable OverhangSet.")
        row, col = self._indexing_helper((idx, side))

        if col == 'V':
            options = [oh for oh in self.valid_v_overhangs[row] if self._not_used_yet(oh)]
            if len(options) == 0:
                self.previous_mutation = None
                return False
            self.previous_mutation = (row, col, self.v_oh_set[row])
            self.v_oh_set[row] = random.choice(options)
            return True

        elif col == 'C':
            options = [oh for oh in self.valid_c_overhangs if self._not_used_yet(oh)]
            if len(options) == 0:
                self.previous_mutation = None
                return False
            self.previous_mutation = (row, col, self.c_oh_set[row])
            self.c_oh_set[row] = random.choice(options)
            return True
        

    def dump_to_dict(self):
        """
        Return dictionary with all critical information about the OverhangSet.
        """

        return {
            'v_seqs': self.v_seqs,
            'c_seq': self.c_seq,
            'v_ids': self.v_ids,
            'c_id': self.c_id,
            'restriction_enzyme': self.restriction_enzyme,
            'order': self.order,
            'is_mutable': self.is_mutable,
            'banned_OHs': self.banned_OHs,
            'v_oh_set': self.v_oh_set.tolist(),
            'c_oh_set': self.c_oh_set.tolist(),
        }


    def save(self, out_path, format):
        """
        Save the OverhangSet to a file in specified format ('csv' or 'json').
        out_path: string, path to output file
        format: string, either 'csv' or 'json'"""

        if format not in ['csv', 'json']:
            raise ValueError("Invalid format specified. Use 'csv' or 'json'.")
        
        if format == 'json':
            with open(out_path, 'w') as f:
                json.dump(self.dump(), f, indent=4)

        elif format == 'csv':
            self.get_OH_set().to_csv(out_path, index=True)
    

    def get_OH_set(self):
        """
        Return a pandas DataFrame with the overhang set, indexed by V_ID.
        Columns are V_ID, V_OH, C_OH (order depends on self.order).
        """

        v_title = self.v_ids[0][:4]
        c_title = self.c_id[:4]

        output = pd.DataFrame()
        output[f"{v_title}_ID"] = self.v_ids
        output[f"{v_title}_OH"] = self.v_oh_set
        output[f"{c_title}_OH"] = self.c_oh_set
        output.set_index(f"{v_title}_ID", inplace=True)

        if self.order == OverhangSet.order_VLeft_CRight:
            output = output[[f"{v_title}_OH", f"{c_title}_OH"]]
        else:
            output = output[[f"{c_title}_OH", f"{v_title}_OH"]]

        return output
    

    # Private helper methods for OverhangSet-------

    def _prepare_valid_overhangs(self):
        """
        Prepare the valid overhang lists for variable and constant regions based on sequences.
        Must be called before random initialization or mutation.
        """

        if not self.is_mutable:
            raise ValueError("Cannot prepare valid overhangs for a non-mutable OverhangSet.")
    
        c_oh_options = OverhangSet._get_OH_from_seq(self.c_seq, self.banned_OHs)
        v_oh_options = [OverhangSet._get_OH_from_seq(v_seq, self.banned_OHs) for v_seq in self.v_seqs]

        #sanity check for constant region overhangs
        all_possible_c_seqs = OverhangSet._dfs_seq_builder(self.c_seq)
        for opt in c_oh_options:
            found = False
            for seq in all_possible_c_seqs:
                if opt in seq:
                    found = True
                    break
            if not found:
                raise ValueError(f"Overhang possibility {opt} not found in constant region sequence.")
            
        #sanity check for variable region overhangs
        for i in range(len(self.v_seqs)):
            all_possible_v_seqs = OverhangSet._dfs_seq_builder(self.v_seqs[i])
            for opt in v_oh_options[i]:
                found = False
                for seq in all_possible_v_seqs:
                    if opt in seq:
                        found = True
                        break
                if not found:
                    raise ValueError(f"Overhang possibility {opt} not found in variable region sequence {self.v_ids[i]}.")

        self.valid_v_overhangs = v_oh_options
        self.valid_c_overhangs = c_oh_options


    def _not_used_yet(self, oh):
        """
        Check if an overhang is not already used in the current set (including reverse complement).
        """
        if oh in self.v_oh_set or oh in self.c_oh_set:
            return False
        if OverhangSet._rev_comp(oh) in self.v_oh_set or OverhangSet._rev_comp(oh) in self.c_oh_set:
            return False
        return True


    def _indexing_helper(self, key):
        """
        Helper function to parse indexing keys for __getitem__ and __setitem__.
        key: tuple of (row, col) where row is int or V_ID string, col is int (0 or 1) or string ('V' or 'C')
        Returns (row_index, col_str) where row_index is int index and col_str is 'V' or 'C'.
        """

        row = key[0]
        col = key[1]

        # if row is a string, find the index using v_ids
        if isinstance(row, str):
            if row not in self.v_ids:
                raise KeyError(f"Variable ID '{row}' not found in OverhangSet.")
            row = self.v_ids.index(row)
        elif isinstance(row, int):
            if row < 0 or row >= len(self):
                raise IndexError(f"Row index {row} out of bounds for OverhangSet with {len(self.v_ids)} entries.")
        else:
            raise TypeError("Row index must be a string (V_ID) or an integer.")

        # if col is a string, determine if it's variable or constant
        v_title = self.v_ids[0][:4]
        c_title = self.c_id[:4]
        if isinstance(col, str):
            if col == "V" or col == v_title:
                col = 'V'
            elif col == "C" or col == c_title:
                col = 'C'
            else:
                raise KeyError(f"Column '{col}' not recognized. Use 'V' or '{v_title}' or 'C' or '{c_title}' if indexing with a string.")
        elif isinstance(col, int):
            if self.order == OverhangSet.order_VLeft_CRight:
                mapper = {0: 'V', 1: 'C'}
            else:
                mapper = {0: 'C', 1: 'V'}

            if col not in [0, 1]:
                raise IndexError(f"Column index {col} out of bounds. Use 0 for left side of vector and 1 for right side of vector.")
            else:
                col = mapper[col] 
        else:
            raise TypeError(f"Column index must be a string ('V' or 'C' or '{v_title}' or '{c_title}') or an integer (0 or 1) for left or right side.")
        
        return row, col
    

    # Dunder methods for OverhangSet-------

    def __eq__(self, other):
        """
        Check equality between two OverhangSet objects.
        Two OverhangSet objects are equal if they have the same v_ids, c_id, restriction_enzyme,
        order, and identical overhang sets.
        """

        if not isinstance(other, OverhangSet):
            return False
        
        if self.v_ids != other.v_ids or self.c_id != other.c_id or self.restriction_enzyme != other.restriction_enzyme:
            return False
        
        if not np.array_equal(self.v_oh_set, other.v_oh_set) or not np.array_equal(self.c_oh_set, other.c_oh_set):
            return False
        
        return True
    

    def __ne__(self, other):
        """
        Check inequality between two OverhangSet objects.
        """
        return not self.__eq__(other)


    def copy(self, is_mutable=None):
        """
        Create a copy of the OverhangSet.
        If is_mutable is specified, the copy will have that mutability. Otherwise, it will match the original.
        """

        mutable = self.is_mutable
        if is_mutable is not None:
            mutable = is_mutable
        new = OverhangSet(self.v_seqs, 
                          self.c_seq, 
                          self.v_ids,
                          self.c_id,
                          self.restriction_enzyme,
                          self.order, 
                          mutable,
                          random_init=False,
                          banned_OHs=self.banned_OHs)
        
        new.v_oh_set = np.copy(self.v_oh_set)
        new.c_oh_set = np.copy(self.c_oh_set)
        return new
    

    def __copy__(self):
        """
        Dunder method for copy.copy().
        """
        return self.copy()
    
    def __str__(self):
        """
        String representation of the OverhangSet.
        """

        v_title = self.v_ids[0][:4]
        c_title = self.c_id[:4]
        output = f"OverhangSet with {len(self)} entries:\n"
        if self.order == OverhangSet.order_VLeft_CRight:
            output = f"{v_title}-{c_title} " + output
        else:
            output = f"{c_title}-{v_title} " + output

        if not self.is_mutable:
            output = "IMMUTABLE " + output
        
        output += f"Restriction Enzyme: {self.restriction_enzyme}\n"
        output += str(self.get_OH_set())
        return output
    

    def __getitem__(self, key):
        row, col = self._indexing_helper(key)
        if col == 'V':
            return self.v_oh_set[row]
        elif col == 'C':
            return self.c_oh_set[row]
        
    
    def __setitem__(self, key, value):
        if not self.is_mutable:
            raise ValueError("Cannot set items in a non-mutable OverhangSet.")
        row, col = self._indexing_helper(key)
        if col == 'V':
            if value not in self.valid_v_overhangs[row]:
                raise ValueError(f"Invalid variable overhang '{value}' for V_ID '{self.v_ids[row]}'.")
            self.v_oh_set[row] = value
        elif col == 'C':
            if value not in self.valid_c_overhangs:
                raise ValueError(f"Invalid constant overhang '{value}' for C_ID '{self.c_id}'.")
            self.c_oh_set[row] = value


    def __len__(self):
        return len(self.v_ids)
    
    
    # Static helper tools to generate valid overhang list-------
    @staticmethod
    def _dfs_seq_builder(ref_aa, primer='', aa_pos=0):
        """
        Depth-first search to build all possible DNA sequences from a reference amino acid sequence.
        ref_aa: string, reference amino acid sequence
        """
        if aa_pos >= len(ref_aa):
            return [primer]
        
        output = []
        for codon_option in OverhangSet._codon_dict[ref_aa[aa_pos]]:
            output.extend(OverhangSet._dfs_seq_builder(ref_aa, primer+codon_option, aa_pos+1))
        return output 


    @staticmethod
    def _unique_4mer_set(seq_list):
        """
        Generate all unique 4-mers for each sequence in the input list seq_list.
        """
        output = set()
        for seq in seq_list:
            for i in range(len(seq) - 3):
                output.add(seq[i:i+4])
        return output

    @staticmethod
    def _rev_comp(seq):
        """
        Reverse complement of a DNA sequence.
        """
        mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join([mapper[nt] for nt in seq[::-1]])

    @staticmethod
    def _comp(seq):
        """
        Complement of a DNA sequence (not reversed).
        """
        mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join([mapper[nt] for nt in seq])

    @staticmethod
    def _is_palindrome(seq):
        """
        Check if a sequence is palindromic (equal to its reverse complement).
        """
        return seq == OverhangSet._rev_comp(seq)

    @staticmethod
    def _remove_revcomp(overhang_set):
        """
        Remove reverse complements from a set of overhangs, keeping only one of each pair.
        """
        new_set = set()
        for overhang in overhang_set:
            if OverhangSet._rev_comp(overhang) not in new_set:
                new_set.add(overhang)
        return new_set

    @staticmethod
    def _get_OH_from_seq(seq, banned=[]):
        """
        Get valid 4bp overhangs from a protein sequence, excluding palindromic and banned overhangs.
        """
        oh_options = OverhangSet._unique_4mer_set(OverhangSet._dfs_seq_builder(seq))
        oh_options = list(filter(lambda x: not OverhangSet._is_palindrome(x), oh_options))
        oh_options = OverhangSet._remove_revcomp(oh_options)
        oh_options = list(filter(lambda x: x not in banned, oh_options))
        return list(oh_options)
    
    _codon_dict = {
            'I': ['ATT','ATC','ATA'],
            'L': ['CTT','CTC','CTA','CTG','TTA','TTG'],
            'V': ['GTT','GTC','GTA','GTG'],
            'F': ['TTT', 'TTC'],
            'M': ['ATG'],
            'C': ['TGT','TGC'],
            'A': ['GCT','GCC','GCA','GCG'],
            'G': ['GGT','GGC','GGA','GGG'],
            'P': ['CCT','CCC','CCA','CCG'],
            'T': ['ACT','ACC','ACA','ACG'],
            'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
            'Y': ['TAT','TAC'],
            'W': ['TGG'],
            'Q': ['CAA','CAG'],
            'N': ['AAT','AAC'],
            'H': ['CAT','CAC'],
            'K': ['AAA','AAG'],
            'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
            'D': ['GAT','GAC'],
            'E': ['GAA','GAG'],
            '*': ['TAA','TAG','TGA']
    }



class ScoringFunction:
    """
    Class to score an OverhangSet based on ligation efficiency data from NEB.
    """

    def __init__(self, ligation_data_dir, weight_coef=0.5):
        """
        ligation_data_dir: directory containing NEB ligation efficiency CSV files for BsaI, BbsI, and BsmBI
        weight_coef: float between 0 and 1 to weight product of efficiencies vs minimum vector efficiency
        """
        if weight_coef < 0 or weight_coef > 1:
            raise ValueError("weight_coef must be between 0 and 1.")
        self.weight_coef = weight_coef
        BbsI = pd.read_csv(os.path.join(ligation_data_dir, 'BbsI_ligation_efficiency.csv')).set_index('Overhang')
        BsaI = pd.read_csv(os.path.join(ligation_data_dir, 'BsaI_ligation_efficiency.csv')).set_index('Overhang')
        BsmBI = pd.read_csv(os.path.join(ligation_data_dir, 'BsmBI_ligation_efficiency.csv')).set_index('Overhang')
        self.ligation_tables = {'BbsI': BbsI, 'BsaI': BsaI, 'BsmBI': BsmBI}
        
    @staticmethod
    def _rev_comp(seq):
        """
        Reverse complement of a DNA sequence.
        """
        mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join([mapper[nt] for nt in seq[::-1]])
    
    def NEB_simple_efficiency_score(self, overhang_list, enzyme):
            """
            Calculate the simple NEB ligation eff score for a list of overhangs and a given enzyme.
            Same as what is calculated on https://ligasefidelity.neb.com
            overhang_list: list of 4bp overhang strings
            enzyme: string, one of 'BsaI', 'BbsI', or 'BsmBI'
            """
            overhang_list_rc = [ScoringFunction._rev_comp(i) for i in overhang_list]
            full_sel = overhang_list + overhang_list_rc
            full_ligation_mat = self.ligation_tables[enzyme].loc[full_sel, full_sel]
            match_ligation_mat = full_ligation_mat.loc[overhang_list, overhang_list_rc]
            return np.prod(np.diag(match_ligation_mat) / np.sum(full_ligation_mat[overhang_list], axis=0))
    
    def per_vector_efficiency(self, oh_set):
        """
        Calculates the probability of successful ligation for each vector (V + C) in the OverhangSet.
        oh_set: OverhangSet object
        """
        all_options = np.unique(np.append(oh_set.v_oh_set, oh_set.c_oh_set))
        all_options_rc = np.array([ScoringFunction._rev_comp(i) for i in all_options])
        enzyme = oh_set.restriction_enzyme
        
        fidelity_table = self.ligation_tables[enzyme].loc[np.append(all_options, all_options_rc), np.append(all_options, all_options_rc)]
        correct_insert_prob = []

        for i in range(len(oh_set)):
            if oh_set.order == OverhangSet.order_VLeft_CRight:
                insert_lh = oh_set.v_oh_set[i]
                vec_rh = oh_set.c_oh_set[i]
            else:
                insert_lh = oh_set.c_oh_set[i]
                vec_rh = oh_set.v_oh_set[i]

            vec_lh = ScoringFunction._rev_comp(insert_lh)
            insert_rh = ScoringFunction._rev_comp(vec_rh)

            l_lig_freq = fidelity_table.loc[vec_lh, insert_lh]
            l_lig_prob = 0.5*(l_lig_freq/np.sum(fidelity_table.loc[vec_lh, :]) + l_lig_freq/np.sum(fidelity_table.loc[:, insert_lh]))
            r_lig_freq = fidelity_table.loc[vec_rh, insert_rh]
            r_lig_prob = 0.5*(r_lig_freq/np.sum(fidelity_table.loc[vec_rh, :]) + r_lig_freq/np.sum(fidelity_table.loc[:, insert_rh]))
            correct_insert_prob.append(l_lig_prob * r_lig_prob)

        return np.array(correct_insert_prob)
    
    def prod_score(self, oh_set):
        """
        Calculate the product of NEB simple efficiency scores for the variable and constant overhang sets.
        Simple efficiency score is calculated separately for V and C overhangs and then multiplied.
        oh_set: OverhangSet object
        """
        return self.NEB_simple_efficiency_score(oh_set.v_oh_set.tolist(), oh_set.restriction_enzyme) * \
               self.NEB_simple_efficiency_score(oh_set.c_oh_set.tolist(), oh_set.restriction_enzyme)
    
    def min_vector_efficiency_score(self, oh_set):
        """
        Return the minimum vector ligation efficiency score across all vectors in the OverhangSet.
        oh_set: OverhangSet object
        """
        return np.min(self.per_vector_efficiency(oh_set))
    
    def TCRAFT_overhang_set_score(self, overhang_set, custom_weight=None):
        """
        Computes the default TCRAFT overhang set score as a weighted combination of the product of
        efficiencies and the minimum vector efficiency score.
        overhang_set: OverhangSet object
        custom_weight: float between 0 and 1 to override the default weight_coef for this calculation
        """
        if custom_weight is None:
            w = self.weight_coef
        else:
            if custom_weight < 0 or custom_weight > 1:
                raise ValueError("Custom weight must be between 0 and 1.")
            w = custom_weight
        return w * self.prod_score(overhang_set) + (1 - w) * self.min_vector_efficiency_score(overhang_set)

    def __call__(self, overhang_set):
        # Dummy example: return random score
        # Replace with real scoring logic (e.g., orthogonality, GC balance, melting temp)
        return self.TCRAFT_overhang_set_score(overhang_set)
    

#---------------MCMC OPTIMIZER------------------------

class MCMC_Optimizer:
    """
    Class to perform MCMC optimization on an OverhangSet using a provided scoring function.
    """

    def __init__(self, oh_set, job_title, score_fxn):
        """
        oh_set: OverhangSet object
        job_title: string, title for the optimization job
        score_fxn: callable, function that takes an OverhangSet and returns a numeric score
        """
        self.job_title = job_title
        self.init_oh_set = oh_set.copy(is_mutable=False)
        self.best_oh_set = oh_set.copy(is_mutable=False)
        self.score_fxn = score_fxn
        self.trace = []
        self.best_score = score_fxn(oh_set)
        
    
    @staticmethod
    def _run_MCMC(oh_set, score_fxn, num_iter, temp, cutoff=1000):
        """
        Run a single MCMC optimization chain.
        oh_set: OverhangSet object (must be mutable)
        score_fxn: callable, function that takes an OverhangSet and returns a numeric score
        num_iter: int, number of MCMC iterations to run
        temp: float, temperature parameter for MCMC
        cutoff: int, number of iterations with no accepted mutations before terminating early
        Returns: best_oh_set, best_score, trace, best_oh_index
        """

        s0 = score_fxn(oh_set)
        trace = [[None, None, None, None, None, s0, temp]]
        no_update = 0
        best_oh_set = oh_set.copy(is_mutable=False)
        best_oh_score = s0
        best_oh_index= 0

        for step in tqdm(range(num_iter)):
            #Step 1: mutate a random position, 30 tries to make a valid mutation
            made_mutation = False
            for i in range(30):
                idx = random.randint(0, len(oh_set)-1)
                side = random.choice(['V', 'C'])
                old_OH = oh_set[idx, side]
                if oh_set.mutate(idx, side):
                    made_mutation = True
                    new_OH = oh_set[idx, side]
                    break
            
            if not made_mutation:
                print('Could not make a valid mutation after 30 tries - overhang set has hit a dead end in optimization, terminating now.')
                break

            #Step 2: Update score and run decision criteria. If s_new > s0, then accept. Otherwise, accept given exp. probability
            # i used log probability to avoid overflow errors
            # update best solution each round
            s_new = score_fxn(oh_set)
            p_acc = (s_new - s0)/temp
            accept_mutation = np.log(random.random()) < p_acc

            if accept_mutation:
                no_update = 0
                s0 = s_new
                if s_new > best_oh_score:
                    best_oh_set = oh_set.copy(is_mutable=False)
                    best_oh_score = s_new
                    best_oh_index = step
            else:
                oh_set[idx, side] = old_OH
                no_update += 1

            trace.append([idx, side, old_OH, new_OH, bool(accept_mutation), s0, temp])

            if no_update >= cutoff:
                print(f'No updates for past {cutoff} iterations, terminating now.')
                break
        
        return best_oh_set, best_oh_score, trace, best_oh_index
    

    def optimize(self, temp, num_iter, random_init, n_seeds, n_processes):
        """
        The main MCMC optimization function. Runs multiple MCMC chains in parallel and updates the best overhang set if a better one is found.
        temp: float, temperature parameter for MCMC
        num_iter: int, number of MCMC iterations per seed
        random_init: bool, whether to randomly initialize each seed or start optimization run from the current best OH set stored in the class
        n_seeds: int, number of parallel MCMC chains to run
        n_processes: int, number of parallel processes to use
        Returns: None (updates internal state of the optimizer)
        """

        print(f"Starting MCMC optimization with random_init={random_init}, {n_seeds} seeds, temp {temp}, and {num_iter} iterations per seed.")

        seeds = [self.best_oh_set.copy(is_mutable=True) for i in range(n_seeds)]
        if random_init:
            for i in seeds:
                i.random_initialize()
            seed_backup = [i.copy(is_mutable=False) for i in seeds]
        
        optim_output = Parallel(n_jobs=n_processes)(delayed(MCMC_Optimizer._run_MCMC)(seed_oh_set, self.score_fxn, num_iter, temp) for seed_oh_set in seeds)
       
        highest_score = optim_output[0][1]
        highest_score_idx = 0
        for i in range(len(optim_output)):
            if optim_output[i][1] > highest_score:
                highest_score = optim_output[i][1]
                highest_score_idx = i

        print(f"Highest score found across {n_seeds} seeds: Seed {highest_score_idx}, {highest_score} at step {optim_output[highest_score_idx][3]}.")

        fig, ax = MCMC_Optimizer._make_optim_plot(optim_output[highest_score_idx][2])
        plt.show()

        if random_init:
            self.best_oh_set = optim_output[highest_score_idx][0]
            self.best_score = optim_output[highest_score_idx][1]
            self.trace = optim_output[highest_score_idx][2][0:optim_output[highest_score_idx][3]]
            self.init_oh_set = seed_backup[highest_score_idx]
            print("Updated best overhang set.")
        else:
            if highest_score > self.best_score:
                self.best_oh_set = optim_output[highest_score_idx][0]
                self.best_score = optim_output[highest_score_idx][1]
                self.trace.extend(optim_output[highest_score_idx][2][1:optim_output[highest_score_idx][3]])
                print("Updated best overhang set.")
            else:
                print("Did not update best overhang set since score did not improve.")

        return
    

    @staticmethod
    def _make_optim_plot(trace):
        """
        Takes an optimization trace (list of lists) and makes a line plot of the score over time colored by temperature.
        trace: list of lists, each inner list is [idx, side, old_OH, new_OH, accepted (bool), score, temp]
        Returns: matplotlib figure and axis objects
        """

        trace_df = pd.DataFrame(trace).iloc[:, 5:]
        trace_df.columns = ['Score', 'Temp']
        trace_df['Step'] = np.arange(0, trace_df.shape[0])
        
        #make a line plot of the score over time colored by temperature, and return this plot (do not show it)
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.lineplot(data=trace_df, x='Step', y='Score', hue='Temp', palette='viridis', ax=ax)
        ax.set_title('MCMC Optimization Trace')
        ax.set_xlabel('Step')
        ax.set_ylabel('Score')
        ax.grid(True)
        plt.tight_layout()
        return fig, ax


    def plot_trace(self, out_path=None):
        """
        Plot the optimization trace stored in the MCMC_Optimizer class instance
        out_path: string or None, if provided the plot will be saved to this path, otherwise it will be shown
        Returns: None
        """

        fig, ax = MCMC_Optimizer._make_optim_plot(self.trace)
        if out_path is not None:
            plt.savefig(out_path)
        else:
            plt.show()
        return

    
    def save(self, out_dir):
        """
        Save the current state of the MCMC_Optimizer to the specified output directory.
        This includes the best overhang set, the optimization trace plot, and a JSON log file.
        out_dir: string, path to output directory
        Returns: None
        """

        try:
            os.mkdir(out_dir)
        except FileExistsError:
            pass
        
        self.best_oh_set.save(os.path.join(out_dir, f"{self.job_title}_Overhang_Set.csv"), format='csv')
        self.plot_trace(os.path.join(out_dir, f"{self.job_title}_MCMC_Score_Trace.png"))

        log_dict = {
            'job_title' : self.job_title,
            'init_oh_set' : self.init_oh_set.dump_to_dict(),
            'best_oh_set' : self.best_oh_set.dump_to_dict(),
            'best_score' : self.best_score,
            'trace': self.trace
        }
        with open(os.path.join(out_dir, f"{self.job_title}_MCMC_Log.json"), 'w') as f:
            json.dump(log_dict, f, indent=4)
        return
    

    @classmethod
    def load_from_json(cls, json_path, ligation_data_dir):
        """
        Load an MCMC_Optimizer instance from a JSON log file.
        json_path: string, path to JSON log file
        ligation_data_dir: string, path to directory containing NEB ligation efficiency CSV files
        Returns: MCMC_Optimizer instance
        """
        
        with open(json_path, 'r') as f:
            log_dict = json.load(f)
        
        init_oh_set = OverhangSet.load_from_dict(log_dict['init_oh_set'])
        best_oh_set = OverhangSet.load_from_dict(log_dict['best_oh_set'])
        score_fxn = ScoringFunction(ligation_data_dir, weight_coef=log_dict.get('weight_coef', 0.5))
        
        optimizer = cls(best_oh_set, log_dict['job_title'], score_fxn)
        optimizer.init_oh_set = init_oh_set
        optimizer.best_score = log_dict['best_score']
        optimizer.trace = log_dict['trace']
        
        return optimizer


if __name__ == "__main__":
    pass