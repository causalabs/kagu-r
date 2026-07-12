test_that("DAG edges are trimmed in physical units", {
  withr::local_pdf(tempfile(fileext = ".pdf"))
  dag <- list(a = character(), b = "a")
  plot <- kagu_plot_dag(dag)

  built <- ggplot2::ggplot_build(plot)
  edge_grob <- plot$layers[[1]]$geom$draw_panel(
    built$data[[1]], built$layout$panel_params[[1]], plot$coordinates
  )
  expect_s3_class(edge_grob, "kagu_dag_edges")
  expect_equal(edge_grob$start_gap, 4.8)
  expect_equal(edge_grob$end_gap, 5.8)
})

test_that("DAG plots without edges still render", {
  withr::local_pdf(tempfile(fileext = ".pdf"))
  plot <- kagu_plot_dag(list(a = character(), b = character()))
  expect_length(plot$layers, 2L)
  expect_silent(ggplot2::ggplotGrob(plot))
})
