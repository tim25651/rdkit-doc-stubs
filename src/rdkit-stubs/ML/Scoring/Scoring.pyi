"""
rdkit.ML.Scoring.Scoring module¶
$Id$
Scoring - Calculate rank statistics
Created by Sereina Riniker, October 2012
after a file from Peter Gedeck, Greg Landrum

param scores: ordered list with descending similarity containingactive/inactive information

param col: column index in scores where active/inactive information is stored
param fractions: list of fractions at which the value shall be calculated
param alpha: exponential weight
"""

def CalcROC(self, scores, col):
    """
    Determines a ROC curve"""
    ...

def CalcAUC(self, scores, col):
    """
    Determines the area under the ROC curve"""
    ...

def CalcRIE(self, scores, col, alpha):
    """
    RIE original definded here:
    Sheridan, R.P., Singh, S.B., Fluder, E.M. & Kearsley, S.K.
    Protocols for Bridging the Peptide to Nonpeptide Gap in Topological Similarity Searches.
    J. Chem. Inf. Comp. Sci. 41, 1395-1406 (2001)."""
    ...

def CalcBEDROC(self, scores, col, alpha):
    """
    BEDROC original defined here:
    Truchon, J. & Bayly, C.I.
    Evaluating Virtual Screening Methods: Good and Bad Metric for the “Early Recognition”
    Problem. J. Chem. Inf. Model. 47, 488-508 (2007).
    ** Arguments**

    scores: 2d list or numpy array0th index representing sample
    scores must be in sorted order with low indexes “better”
    scores[sample_id] = vector of sample data

    col: intIndex of sample data which reflects true label of a sample
    scores[sample_id][col] = True iff that sample is active

    alpha: floathyper parameter from the initial paper for how much to enrich the top

    Returnsfloat BedROC score"""
    ...

def CalcEnrichment(self, scores, col, fractions):
    """
    Determines the enrichment factor for a set of fractions"""
    ...

def CalcRIE(self, scores, col, alpha):
    """
    RIE original definded here:
    Sheridan, R.P., Singh, S.B., Fluder, E.M. & Kearsley, S.K.
    Protocols for Bridging the Peptide to Nonpeptide Gap in Topological Similarity Searches.
    J. Chem. Inf. Comp. Sci. 41, 1395-1406 (2001)."""
    ...

def CalcROC(self, scores, col):
    """
    Determines a ROC curve"""
    ...
