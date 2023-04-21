# DELTA

DELTA (DNA Encoded Library Topological Assignment) is a system for classificaion of DNA encoded libraries (DELs) according to the topological arrangement of thier building blocks (BBs). This system assigns topological classification to a library using the following algorithm
The following repository is included as supporting information for **LINK**
DELTA (DNA Encoded Library Topological Assignment) is a system for classification of DNA encoded libraries (DELs) according to the topological arrangement of their building blocks (BBs). This system assigns topological classification to a library as single string using the following algorithm:
1)	Identify the longest path to the DNA attachment point.
2)	Moving along the main chain, alphabetically assign each structural element in ascending order. Adjacent elements share a connection it is a branching element.
3)	If a branch point in the main chain is reached, the branching element is included after the branch point in parentheses (X) using next sequential letter before continuing on with the main chain.
4)	Cycles are handled by defining a cyclic closure between two elements using two exclamation points. From left to right, the first cycle point is considered the cycle start and indicated by placing the exclamation point before the element (!X). The next cycle point is considered the cycle end and is indicated by putting the exclamation point after the element (X!).
