# Version 1.0.0.10
## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

* Fix the problem with sprintf on M1mac (https://www.stats.ox.ac.uk/pub/bdr/sprintf.txt)


# Version 1.0.0.9
## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* extended package description
* removed Rd files for non-export functions
* clean up the author definitions
* rephrase package description

## Reivewer comment
Uwe Ligges:
```
Thanks, we see:



   The Description field should not start with the package name,
     'This package' or similar.

Please also single quote software names such as 'rSpectral' in the
Description field
```
### Submission comments:

Address all comments, in particular add single quote to package name in the
DESCRIPTION file.

# Version 1.0.0.8
### Submission comments:

Please write references in the description of the DESCRIPTION file in
the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: authors (year) <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Please always add all authors, contributors and copyright holders in the
Authors@R field with the appropriate roles.
 From CRAN policies you agreed to:
"The ownership of copyright and intellectual property rights of all
components of the package must be clear and unambiguous (including from
the authors specification in the DESCRIPTION file). Where code is copied
(or derived) from the work of others (including from R itself), care
must be taken that any copyright/license statements are preserved and
authorship is not misrepresented.
Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
authors of such code. Alternatively, the ‘Author’ field should list
these authors as contributors.
Where copyrights are held by an entity other than the package authors,
this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’
field, or using a ‘Copyright’ field (if necessary referring to an
inst/COPYRIGHTS file)."
e.g.: "Written by Mark Newman 11 AUG 06 "
Please explain in the submission comments what you did about this issue.

#### Responce
1. I've updated references in Description according to the format.
2. All authors and contributors are listed in the Description file with 
appropriate comments desrcibing their roles.
