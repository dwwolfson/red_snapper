Readme for parameterized reports

This is an explanation of the process for varying grid cell size with the ctmcmove method of Hanks et al. 2015.

A single file, 'ctmc_parameterized.Rmd', can do serve several functions. By running a simple function you can specify the size you'd like to use and the file will create a html showing grouped ctmc model results as well as individual coefficients, and also save a version of the data just prior to fitting the glm model.

As of 4/11/2021, we're looking at just 3 snappers for a comparision example with iSSFs (num 9,12,16). I made htmls for 5-50 m^2 by 5. My next step is to import all data intermediate files, fit each set, and plot the coefficients varying by cell size.

In order to run this .Rmd file using parameter specification first run the last chunk:  
```  
  render_report<-function(resolution, title){
    rmarkdown::render(
      here::here("scripts/parameterized_reports/ctmc_parameterized.Rmd"), params=list(
        resolution=resolution,
        title=title
      ),
      output_file = paste0("CTMC 3 snapper example grouped and individual summary for ", resolution, " m^2 grid cell   size.html")
      )
  }  
```  
  
  Then if you want to run the script at 10m resolution, just run render_report(10, "add title here") in the console.