#11/3/2025
#Aland Chan
#This document lists R functions written over the years
#Mainly to process GIS data
#Upload to github

#1. Housekeeping GIS ####

#1.1 Function to mosaic buffered tiles 

# Currently, the most efficient way to mosaic tiles is to use the gdalbuildvrt function followed by gdalwarp. However, a problem with this workflow (and GDAL in general) is that it always transfers values from the first raster to the output when two rasters overlap. This is unsatisfactory as it defeats the purpose of avoiding  

gdal_buffered_mosaic<-function(input_tilespath=NA,output_filename,tilelist=NA,buffer=NA,ncores,rm_tmp=T){
  
  require(gdalUtilities)
  
  #If there is a buffer to remove
  
  if(is.na(buffer)!=T){
    
    require(foreach)
    require(doParallel)
    
    if(is.na(tilelist)[1]==T){
      
      #We specified the filepath to the tiles not the tilelist itself
      #First list the tif files
      
      tilelist<-list.files(input_tilespath,full.names = T,all.files = F,recursive = T)
      tilelist<-grep(".tif",tilelist,value = T)
      
    }
    
    #Create a directory to store the temporary cropped images
    
    dir.create(paste0(dirname(output_filename),"/temp_cropped_dir"))
    temp_cropped_dir<-paste0(dirname(output_filename),"/temp_cropped_dir")
    
    #Find the corners of the input rasters
    
    c2<-makeCluster(ncores)
    registerDoParallel(c2)
    
    new_ullr_list<-foreach(x=tilelist) %dopar% {
      
      info<-system(paste0("gdalinfo ",x),intern = T)
      
      ul<-grep("Upper Left",info,value = T)
      lr<-grep("Lower Right",info,value = T)
      
      ulx<-regmatches(ul,regexpr("\\([[:space:]]*[0-9.]+",ul))
      ulx<-as.numeric(regmatches(ulx,regexpr("[0-9.]+",ulx)))
      ulx<-ulx+buffer/2
      
      
      uly<-regmatches(ul,regexpr("[0-9.]+[[:space:]]*\\)",ul))
      uly<-as.numeric(regmatches(uly,regexpr("[0-9.]+",uly)))
      uly<-uly-buffer/2
      
      lrx<-regmatches(lr,regexpr("\\([[:space:]]*[0-9.]+",lr))
      lrx<-as.numeric(regmatches(lrx,regexpr("[0-9.]+",lrx)))
      lrx<-lrx-buffer/2
      
      lry<-regmatches(lr,regexpr("[0-9.]+[[:space:]]*\\)",lr))
      lry<-as.numeric(regmatches(lry,regexpr("[0-9.]+",lry)))
      lry<-lry+buffer/2
      
      corners<-c(ulx,uly,lrx,lry)
      
      return(corners)
      
    }
    
    stopCluster(c2)
    
    #Check extent, if the buffer is larger than the raster, skip the file
    
    check<-sapply(new_ullr_list,function(x){
      return((x[3]-x[1])>0 & (x[2]-x[4])>0)
    })
    
    r_list_checked<-c(1:length(tilelist))[which(check==T)]
    
    #Crop the remaining rasters
    
    c2<-makeCluster(ncores)
    registerDoParallel(c2)
    
    foreach(i=r_list_checked,.packages="gdalUtilities") %dopar% {
      
      gdal_translate(
        src_dataset = tilelist[i],
        dst_dataset = paste0(temp_cropped_dir,"/",gsub(".tif","_cropped.tif",basename(tilelist[i]))),
        projwin = new_ullr_list[[i]]
      )
      
    }
    
    stopCluster(c2)
    
    #List the cropped rasters
    
    cr_fl<-list.files(temp_cropped_dir,full.names = T,all.files = F,recursive = T)
    
    #Mosaic the cropped tiles
    
    gdalbuildvrt(gdalfile = cr_fl,
                 output.vrt = paste0(temp_cropped_dir,"/temp_mosaic.vrt"))
    
    gdalwarp(srcfile = paste0(temp_cropped_dir,"/temp_mosaic.vrt"),
             dstfile = output_filename,
             co = "COMPRESS=NONE",
             wo = "NUM_THREADS=val/ALL_CPUS",
             multi = T)
    
    #Remove temp files
    
    if(rm_tmp==T){
      
      unlink(temp_cropped_dir,recursive = T)
      
    }
    
  } else {
    
    #No buffer to remove
    
    tilelist<-list.files(input_tilespath,full.names = T,all.files = F,recursive = T)
    tilelist<-grep(".tif",tilelist,value = T)
    
    gdalbuildvrt(gdalfile = tilelist,
                 output.vrt = paste0(dirname(output_filename),"/temp_mosaic.vrt"))
    
    gdalwarp(srcfile = paste0(dirname(output_filename),"/temp_mosaic.vrt"),
             dstfile = output_filename,
             co = c("TILED=YES","COMPRESS=NONE"),
             wo = "NUM_THREADS=val/ALL_CPUS",
             multi = T)
    
    if(rm_tmp==T){
      
      unlink(paste0(dirname(output_filename),"/temp_mosaic.vrt"))
      
    }
    
  }
  
  
}

# 1.2 Split tiles 

raster_split<-function(raster_fp,tiles_factor,tiles_dir,ncores){
  
  library(foreach)
  library(doParallel)
  library(gdalUtilities)
  
  #Read the first raster
  
  raster<-raster::stack(raster_fp)
  
  #Create index
  
  nrows<-raster::nrow(raster)
  ncols<-raster::ncol(raster)
  
  row_width<-round(nrows/tiles_factor)
  col_width<-round(ncols/tiles_factor)
  
  x_seq<-seq(from=0,to=ncols,by=col_width)
  y_seq<-seq(from=0,to=nrows,by=row_width)
  
  rm(raster)
  
  c2<-makeCluster(ncores)
  registerDoParallel(c2)
  
  foreach(x=x_seq,.inorder = F,.packages = "gdalUtilities") %dopar% {
    
    if (x+col_width <= ncols){
      x_width<-col_width
    } else {
      x_width<-ncols-x
    }
    
    #Now loop through all the y values for this x
    
    for (y in y_seq){
      
      if (y+row_width <= nrows){
        y_width<-row_width
      } else {
        y_width<-nrows-y
      }
      
      fp<-paste0(
        tiles_dir,"/",
        gsub(".tif$",
             paste0("_",x,"_",y,".tif"),
             basename(raster_fp))
      )
      
      gdal_translate(src_dataset = raster_fp,
                     dst_dataset = fp,
                     srcwin = c(x,y,x_width,y_width))
      
    }#End of y loop
    
  } #End of foreach loop
  
  stopCluster(c2)
}



# 1.9 gdaltilewarp

gdaltilewarp<-function(tiles_factor,buffer,ncores,rm_tmp=T,srcfile,dstfile,...){
  
  require(foreach)
  require(doParallel)
  require(raster)
  require(gdalUtilities)
  
  #Read the first raster
  
  raster<-raster::stack(srcfile)
  
  #First split the raster into tiles_factor*ncores number of tiles
  
  #Create index
  
  ext<-extent(raster)
  ext<-c(ext@xmin,ext@xmax,ext@ymin,ext@ymax)
  
  nrows<-as.numeric(as.character(ext[4]))-as.numeric(as.character(ext[3]))
  ncols<-as.numeric(as.character(ext[2]))-as.numeric(as.character(ext[1]))
  
  row_width<-round(nrows/tiles_factor)
  col_width<-round(ncols/tiles_factor)
  
  x_seq<-seq(from=as.numeric(as.character(ext[1])),to=as.numeric(as.character(ext[2])),by=col_width)
  y_seq<-seq(from=as.numeric(as.character(ext[3])),to=as.numeric(as.character(ext[4])),by=row_width)
  
  rm(raster)
  
  #Create temporary directory for tiles
  
  dir.create(paste0(dirname(dstfile),"/temp_warp_tiles"))
  dir.create(paste0(dirname(dstfile),"/temp_warp_tiles_output"))
  
  c2<-makeCluster(ncores)
  registerDoParallel(c2)
  
  foreach(x=x_seq,.packages = "gdalUtilities") %dopar% {
    
    if (x+col_width+buffer <= as.numeric(as.character(ext[2]))){
      end_x<-x+col_width+buffer
    } else {
      end_x<-as.numeric(as.character(ext[2]))
    }
    
    #Now loop through all the y values for this x
    
    for (y in y_seq){
      
      if (y+row_width+buffer <= as.numeric(as.character(ext[4]))){
        end_y<-y+row_width+buffer
      } else {
        end_y<-as.numeric(as.character(ext[4]))
      }
      
      tile_fp<-paste0(dirname(dstfile),"/temp_warp_tiles/"
                      ,"tile_",x,"_",y,".tif")
      
      #Create the tile
      
      gdal_translate(src_dataset = srcfile,
                     dst_dataset = tile_fp,
                     r = "bilinear",
                     projwin = c(x,end_y,end_x,y))
      
      #Warp the tile
      # 
      # gdalwarp(srcfile = tile_fp,
      #    dstfile = paste0(dirname(dstfile),
      #                     "/temp_warp_tiles_output/",
      #                     basename(tile_fp)),
      #    t_srs = "EPSG:2326")
      
      gdalwarp(srcfile = tile_fp,
               dstfile = paste0(dirname(dstfile),
                                "/temp_warp_tiles_output/",
                                basename(tile_fp)),
               ...)
      
      
    } #y for loop end
    
  }#x foreach loop end
  
  stopCluster(c2)
  
  #Mosaic the tiles back together
  
  gdal_buffered_mosaic(input_tilespath=paste0(dirname(dstfile),"/temp_warp_tiles_output"),
                       output_filename=dstfile,
                       buffer=buffer,
                       ncores=ncores,
                       rm_tmp=T)
  
  #Remove temp files
  
  if(rm_tmp==T){
    
    unlink(paste0(dirname(dstfile),"/temp_warp_tiles"),recursive = T)
    unlink(paste0(dirname(dstfile),"/temp_warp_tiles_output"))
    unlink(paste0(dirname(dstfile),"/gdaltileswarp.vrt"))
    
  }
  
} #End of gdaltileswarp function



# 1.14 st_find_intersect 

st_find_intersect<-function(sf1,sf2,intersect=T){
  
  require(sf)
  
  intersects_list<-st_intersects(sf1,sf2)
  
  int_vec<-sapply(intersects_list,length) 
  
  #int_vec is a vector showing how many polygons in sf2 a polygon in sf1 intersects with
  #If 0, no intersection, if >0, intersection
  
  if(intersect==T){
    
    sf1_intersecting<-sf1[which(int_vec>0),]
    return(sf1_intersecting)
    
  } else {
    
    sf1_not_intersecting<-sf1[which(int_vec==0),]
    return(sf1_not_intersecting)
    
  }
  
  
}#End of st_find_intersect function

# 1.16 MatchExtRes 


MatchExtRes<-function(source_r_fp,target_r_fp,output_r_fp,resampling="bilinear",thickness=NA,terra=F,srcnodata="None",a_srs = NA){
  
  
  require(gdalUtilities)
  
  #Set the environmental variables to allow gdal to run in R
  
  # Sys.setenv(PROJ_LIB="C:/Program Files/GDAL/projlib")
  
  #Get number of layers if not specified
  
  if(terra == T){
    
    #Terra method to create new raster
    
    require(terra)
    
    src_r<-terra::rast(source_r_fp)
    src_nlyr<-nlyr(src_r)
    rm(src_r)
    
    tr<-terra::rast(target_r_fp)
    
    er<-terra::rast(ext(tr),
                    resolution = res(tr),
                    crs = crs(tr),
                    nlyrs = src_nlyr)
    
    terra::writeRaster(er,output_r_fp)
    
  } else {
    
    require(raster)
    
    #gdal_create method to create new raster
    
    if(is.na(thickness)==T){
      
      temp_raster<-raster::stack(source_r_fp)
      thickness<-raster::nlayers(temp_raster)
      rm(temp_raster)
      
    }
    
    #Create a string in the format "-b 1 -b 2 -b 3 " to input into gdal_translate
    
    temp_string<-1:thickness
    temp_string<-paste0("-b ",temp_string," ")
    temp_string<-paste0(temp_string,collapse="") # "-b 1 -b 2 -b 3 -b 4 "
    
    #gdal_create to create a new file
    #Same extent and resolution as the target with "thickness" number of bands
    #Rename it as the output filename
    
    if(is.na(a_srs)==T){
      
      system(paste0(
        "gdal_create ",
        "-if ",target_r_fp," ",
        "-bands ",thickness," ",
        output_r_fp
      ))
      
    } else {
      
      system(paste0(
        "gdal_create ",
        "-a_srs ",a_srs," ",
        "-if ",target_r_fp," ",
        "-bands ",thickness," ",
        output_r_fp
      ))
      
    }
    
  }
  
  # gdal_translate(src_dataset = target_r_fp,
  #                dst_dataset = output_r_fp,
  #                b = 1:thickness)
  
  # system(paste0("gdal_translate ",
  #               temp_string,
  #               target_r_fp," ",
  #               output_r_fp))
  
  #Warp the source data onto the newly created file
  
  gdalwarp(srcfile = source_r_fp,
           dstfile = output_r_fp,
           r = resampling)
  
  # system(paste0("gdalwarp ",
  #               "-r ",resampling," ",
  #               "-srcnodata None ",
  #               source_r_fp," ",
  #               output_r_fp))
  
}



# 1.17 Visualise binned means

#This offers a quick way to visualise the relationship y~x without resorting to complex density plots.

binned_plot<-function(x_var,y_var,n_bins = 12){
  
  require(dplyr)
  require(ggplot2)
  
  #Create dataframe
  
  df<-data.frame(x = x_var,
                 y = y_var)
  
  #Define bins using quantiles
  
  cut_vec<-seq(from = 0, to = 1, by = 1/n_bins)
  cut_vec<-quantile(x_var,probs = cut_vec)
  
  #Cut the x_var into groups
  
  df$x_groups <- cut(df$x,breaks = cut_vec,include.lowest = T)
  
  #Summarise
  
  df<-df %>% 
    group_by(x_groups) %>% 
    summarise(binned_x = mean(x), mean_y = mean(y)) %>% 
    ungroup()
  
  #Plot
  
  p<-ggplot(df,aes(x=binned_x,y=mean_y)) + 
    geom_point()
  
  return(p)
  
}


### 1.18 This function renames geometry column of sf

rename_geometry <- function(g, name){
  current = attr(g, "sf_column")
  names(g)[names(g)==current] = name
  st_geometry(g)=name
  g
}


#This function gets ullr from a filepath to a raster
#We use gdal to avoid loading things into R

get_ullr<-function(raster_fp){
  
  require(gdalUtilities)
  
  gdalinfo_lines<-gdalinfo(raster_fp,quiet = T)
  
  # Split the output into lines
  gdalinfo_lines <- strsplit(gdalinfo_lines, "\n")[[1]]
  
  gdalinfo_lines<-grep("Upper Left|Lower Right",
                       gdalinfo_lines,
                       value = T)
  
  
  ul<-grep("Upper Left",
           gdalinfo_lines,
           value = T)
  
  ul_xy<-regmatches(ul, gregexpr("\\(.*?\\)", ul))[[1]][1]
  ul_xy<-gsub("[( )]","",ul_xy)
  ul_xy<-strsplit(ul_xy,split = ",")
  ul_x<-ul_xy[[1]][1]
  ul_y<-ul_xy[[1]][2]
  
  lr<-grep("Lower Right",
           gdalinfo_lines,
           value = T)
  
  lr_xy<-regmatches(lr, gregexpr("\\(.*?\\)", lr))[[1]][1]
  lr_xy<-gsub("[( )]","",lr_xy)
  lr_xy<-strsplit(lr_xy,split = ",")
  lr_x<-lr_xy[[1]][1]
  lr_y<-lr_xy[[1]][2]
  
  ullr<-c(ul_x,ul_y,lr_x,lr_y)
  ullr<-as.numeric(ullr)
  return(ullr)
} #End of get_ullr function




#Function to trim df and assign bin index
#Takes a dataframe, column for binning, number of pixels per bin
#Given a certain bin size (n), it is likely that nrow(df) is not divisible by n
#So we need to trim the df and remove some rows to ensure that each bin contains n pixels
#We trim equally from the smallest and largest values

df_bin_index<-function(df,column,n){
  
  #Order by the column selected
  
  df<-df[order(column),]
  
  #Calculate remainder and trim df
  
  remainder<-nrow(df) %% n
  
  if(remainder==1){
    df<-df[-1,]
  } else if (remainder>1){
    rows_to_remove<-c(1:floor(remainder/2),
                      (nrow(df)-ceiling(remainder/2)+1):nrow(df))
    df<-df[-rows_to_remove,]
  }
  
  #Now we could assign groups
  
  df$tgroup<-rep(1:(nrow(df)/n),each = n)
  
  return(df)
  
} #End of df_bin_index function



#Alternative function for st_intersects
#Returns a dataframe with four columns:
#1. Row name of sf1 feature
#2. Row name of sf2 feature that overlaps with the sf1 feature
#3. Percentage of sf1 feature area that is overlapping
#4. Percentage of sf2 feature area that is overlapping

st_area_intersects<-function(sf1,sf2){

  require(sf)
  
  intersections<-st_intersects(sf1,sf2)
  
  intersections<-lapply(1:length(intersections),function(sf1_feature_num){
    
    row_num<-intersections[[sf1_feature_num]]
    
    if(length(row_num)==0){
      
      tdf<-data.frame(
        sf1_row = sf1_feature_num,
        sf2_ovl_row = NA,
        sf1_ovl_perc = NA,
        sf2_ovl_perc = NA
      )
      
      return(tdf)
      
    } else {
      tsf1<-sf1[sf1_feature_num,]
      tsf2<-sf2[row_num,]
      t_ovlp<-suppressWarnings(st_intersection(tsf1,tsf2))
      ovl_area<-as.numeric(st_area(t_ovlp))
      
      tdf<-data.frame(
        sf1_row = sf1_feature_num,
        sf2_ovl_row = row_num,
        sf1_ovl_perc = ovl_area/as.numeric(st_area(tsf1)),
        sf2_ovl_perc = ovl_area/as.numeric(st_area(tsf2))
      )
      
      return(tdf)
    }
    
  }) #End of lapply
  
  intersections<-do.call(rbind,intersections)
  }


#Function to create file path dataframes
fpdf_create<-function(filepath,
                      grep_pattern,
                      id_pattern,
                      id_name = "id",
                      full.name = T,
                      recursive = F){

  require(stringr)

  file_list<-list.files(filepath, full.name = full.name, recursive = recursive)
  file_list<-grep(grep_pattern,file_list,value = T)
  
  file_list<-data.frame(
    filepath = file_list,
    attribute = str_extract(file_list,id_pattern)
  )

  names(file_list)[2]<-id_name

  return(file_list)
}



# Function to create bounding boxes from raster files.
# The advantage of this function is that rasters do not need to be read into memory
# The input is a vector of filepaths, the output is an sf object

raster_bbox_to_sf<-function(filepaths,crs,id_pattern="none"){

  require(sf)
  require(dplyr)
  require(gdalUtilities)
  require(stringr)
  
  bbox_list<-lapply(filepaths,function(fp){

    ullr<-get_ullr(fp)
    ul<-c(ullr[1],ullr[2])
    lr<-c(ullr[3],ullr[4])
    
    # Create a rectangular polygon (order matters: counterclockwise)
    rect_coords <- matrix(
      c(ul[1], ul[2],  # Top-left
        lr[1], ul[2],  # Top-right
        lr[1], lr[2],  # Bottom-right
        ul[1], lr[2],  # Bottom-left
        ul[1], ul[2]), # Close the polygon (back to top-left)
      ncol = 2, byrow = TRUE
    )
    
    # Convert to an sf polygon
    rectangle <- st_sfc(st_polygon(list(rect_coords)), crs = crs)  # Specify CRS (WGS84)
    
    # Create an sf object
    rectangle_sf <- st_sf(geometry = rectangle)
    
    rectangle_sf<-rectangle_sf %>% 
      mutate(filepath = fp,.before = 1)

    if(id_pattern!="none"){

      retangle_sf<-retangle_sf %>%
        mutate(id = str_extract(fp,id_pattern))
      
    }
    
    return(rectangle_sf)
    
  }) #End of lapply

  #Bind the list of sfs back together

  bbox_list<-do.call(rbind,bbox_list)

  return(bbox_list)
  
}
