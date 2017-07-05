context("annotateTraces")

annotatedTraces <- annotateTraces(traces=examplePeptideTracesUnannotated,
                                  trace_annotation=exampleTraceAnnotation,
                                  traces_id_column = "protein_id",
                                  trace_annotation_id_column = "Entry",
                                  trace_annotation_mass_column = "Mass",
                                  uniprot_mass_format = TRUE,
                                  replace_whitespace = TRUE)

test_that("Output format",{
  # test if untouched traces parts stay identical
  testthat::expect_identical(examplePeptideTracesUnannotated$traces, annotatedTraces$traces)
  testthat::expect_identical(examplePeptideTracesUnannotated$trace_type, annotatedTraces$trace_type)
  testthat::expect_identical(examplePeptideTracesUnannotated$fraction_annotation, annotatedTraces$fraction_annotation)
  testthat::expect_identical(examplePeptideTracesUnannotated$trace_annotation$id, annotatedTraces$trace_annotation$id)
  # test if annotation was successful
  testthat::expect_true(all(names(exampleTraceAnnotation) %in% c(names(annotatedTraces$trace_annotation),"Entry")))
  })
