# assign

This program assigns subjects to professors, taking into consideration their preferences. The program uses simulated annealing to optimise a cost function associated with an assignment.

## The program takes two inout files.

### Preferences file (e.g., prefs.txt)

Each line of the preferences file has: name, min_hours, max_hours, 6 preferences, where min_hours (max_hours) is the minimum (maximum) number of hours per week the professor will be assigned. Each preference consists of a weignt (positive integer) and a combination of at most 3 subjects the professor would like to teach. Lower weight for a preference signifies higher preference for the combination of subjects. Sum of the weights of the 6 preferences is to be 15.

### Disciplinas data (e.g., disc.txt)

Each line has 4 fields: id (subject id), code (code of the turma/section for the hours at which a section is taught), num_sec (number of sections with the given subject id and code), credits (the number of credit hours per week). In a single preference, two subjects cannot have the same code, since the code signifies the hours at which a turma (section) is taught,

## To run the program:

assign -p prefs.txt -d disc.txt -z <number>

The parameter passed with -z option is a seed for the pseudo random number generator. This parameter is optional, but it helps if we want to repeat the results.

In principle, the program can be adapted for many kinds of allocation problems, but at the moment it is not generic enough. This will be improved in the future versions.
