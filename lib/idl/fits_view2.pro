function make_contour_levels, range, nlevels,  log

    i = lindgen(nlevels)

    if log eq 0 then begin  ; linear levels
        dz = (range[1] - range[0]) / nlevels
        levels = range[0] + i*dz
    end else begin          ; log levels
        di = alog10(range[1]/range[0])/(nlevels-1)
        levels = range[0] * 10D^(i*di)
    end

    print, "Made "+strn(nlevels)+" contours from "$
        +strn(range[0])+" to "+strn(range[1])+" : "

    return, levels
end

function arrayfy, var, N

	if N_elements(var) eq 1 then $
		var = replicate(var, N)

	return, var
end

pro set_txtct, txtct, i, img_res, img

	if txtct[i] ne -1 then $ 
		return 
            	
	avg_size = [0.1, 0.25]*img_res  ; [pix] area on img to get color from
   	avg = mean(img[0:avg_size[0], img_res-avg_size[0]:*])
            
	if avg gt 128 then $
		txtct[i] = 0      $
	else $
		txtct[i] = 255

	return 
end        

; plot 2D vectorfield whithout arrows
pro vecplot, u, v, Length = length, _REF_EXTRA=extra

        on_error,2       ;Return to caller if an error occurs

        if not keyword_set(title)   then title = ''
        if not keyword_set(length)  then length = 1.

        s = size(u)
        t = size(v)
        
        x = findgen(s[1]) 
        y = findgen(s[2])  

        mag = sqrt(u^2 + v^2)             
                
        good = where(mag lt 1e30, nbad) ; cutoff

        ugood = u[good]
        vgood = v[good]
 
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)

        x_step=(x1-x0)/(s[1]-1.)   
        y_step=(y1-y0)/(s[2]-1.)   

        maxmag = max([max(abs(ugood/x_step)) $
                     ,max(abs(vgood/y_step))])
       
        sin_a = length * (ugood/maxmag)
        cos_a = length * (vgood/maxmag)

        x_b0 = x0-x_step
        x_b1 = x1+x_step
        y_b0 = y0-y_step
        y_b1 = y1+y_step
        
        plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            _EXTRA = extra

      
        for i=0L,n_elements(good)-1 do begin     ;Each point

                x0 = x[good[i] mod s[1]]   ;get coords of start & end
                dx = sin_a[i]
                x1 = x0 + dx
                x0 -= 0.5*dx    ; place vector in the middle of pixel   
                
                y0 = y[good[i] / s[1]]
                dy = cos_a[i]
                y0 -= 0.5*dy
                y1 = y0 + dy
                
                plots,[x0,x1,x1,x1,x1], [y0,y1,y1,y1,y1], $
                      _EXTRA=extra
        end
end

; always extend vertically
function find_layout, Nimg
    
    XYimg = make_array(2,/float)

    XYimg[0] = ceil(sqrt(Nimg))
    XYimg[1] = ceil(Nimg/XYimg[0])

    return, XYimg
end

; find+set range
function set_range, img, range=range, silent=silent

    if range[0] eq range[1] then begin
            range = minmax(img)

            if not keyword_set(silent) then $
                print, "Setting auto range : ", range
    end             
   
    bad = where(img lt range[0],cnt)
    if cnt ne 0 then $
        img[bad] = range[0] 
    
    bad = where(img gt range[1],cnt)
    if cnt ne 0 then $
        img[bad] = range[1]    

    return, img
end

; apply log & scale
function apply_scale, img, range, log, silent=silent

    if log ne 0 then begin 
        img = alog10(img)
        range = alog10(range)
    end

    if not keyword_set(silent) then $
        print, 'Scale : ', range, log

    img = 255*(img[*,*]-range[0])/(range[1]-range[0])

     return, img
end

function zoom_image, img, zoom, silent=silent
 
    res = (size(img))[1]

    idx_start = res*(1/2. - 1./(2*zoom) )
    idx_stop = res*(1/2. + 1./(2*zoom) )-1
    
    if not keyword_set(silent) then $
        print, 'Zooming factor ', zoom

    return, img[idx_start:idx_stop,idx_start:idx_stop]
end

; plot a series of img. all sizes in cm
; internally all input vars are per img, ie arrays
pro fits_view2, fin, range=range, ext=ext, fout=fout, unit=unit, $
    colmap=colmap, log=log, barthick=barthick,vec_plot=vec_plot, $
    nocolbar=nocolbar, factor=factor, pdf=pdf, text=text, frame_plot=frame,  $
	frame_charsize=frame_charsize, smooth_fwhm=smooth_fwhm,silent=silent, $
	std_colmap=std_colmap, zoom=zoom, bartext=bartext, xythick=xythick, $ 
    invert_colmap=invert_colmap, helix_values=helix_values, xyimg=xyimg, $
    txtct=txtct, bar_nticks=bar_nticks, frame_nminor=frame_nminor, cube=cube, $ 
	$ ; contour over plot 
	cont_fin=cont_fin, nlevels=nlevels,cont_smooth=cont_smooth, $
    cont_thick=cont_thick, cont_factor=cont_factor, cont_range=cont_range, $
    cont_label=cont_label, cont_levels=cont_levels, cont_colmap=cont_colmap, $
    cont_charsize=cont_charsize, cont_ext=cont_ext, cont_log=cont_log, $
	cont_inv=cont_inv 

	;on_error, 2

    if not keyword_set(fin) then begin
        print, 'Display multiple images with one colorbar'
        print, 'Usage: '
        print, 'fits_view2, fin, range=range, ext=ext, fout=fout,'
        print, '            unit=unit, colmap=colmap, log=log, std_colmap= '
        print, '            barthick=barthick,vec_plot=vec_plot, factor=factor,'
        print, '            pdf=pdf, text=text, smooth_fwhm=smooth_fwhm'
        print, '            silent=silent,nocolbar=nocolbar, con_fin=cont_fin'
        return
    end

    if not keyword_set(ext)         then ext = 0
    if not keyword_set(fout)        then fout = './fview2.eps'
    if not keyword_set(colmap)      then colmap = -1
    if not keyword_set(log)         then log = 0
    if not keyword_set(barthick)    then barthick = 15
    if not keyword_set(factor)      then factor = 1
    if not keyword_set(text)        then text = ' '
    if not keyword_set(img_res)     then img_res = 1024 ; [pixels], good compromise
    if not keyword_set(range)       then range = double([[0],[0]])
    if not keyword_set(frame)       then frame = 0  
    if not keyword_set(std_colmap)  then std_colmap = 0
    if not keyword_set(zoom)        then zoom = 1
    if not keyword_set(xythick)     then xythick = 3
    if not keyword_set(cont_ext)    then cont_ext = 0
    if not keyword_set(cont_log)    then cont_log = 0
    if not keyword_set(cont_label)  then cont_label = 0
    if not keyword_set(cont_factor) then cont_factor = 1
    if not keyword_set(cont_thick)  then cont_thick = 1
    if not keyword_set(cont_range)  then cont_range = double([[0],[0]])
    if not keyword_set(cont_colmap) then cont_colmap = -1
    if not keyword_set(cont_inv) 	then cont_inv = 0
    if not keyword_set(txtct) 		then txtct = -1

    if not keyword_set(silent) then begin
        print, 'Fits_view2 - IDL fits plotting'
		print, 'Input from ', fin
        print, 'O	utput to ',fout
    end
    
    ; set parameters
    img_res = 1024              ; [pixels], qualitiy & filesize compromise

    Nimg = n_elements(fin)      ; number of img to be shown
    
    if Nimg gt 1 then begin ; make single parameters into arrays, if necessary

        if N_elements(range) eq 2 then begin
            range = [[replicate(range[0], Nimg)],$
                     [replicate(range[1], Nimg)] ] 
            print, "WARNING: Autorange with multiple images !!"
        end

		if N_elements(cont_range) eq 2 then begin
            cont_range = [[replicate(cont_range[0], Nimg)],$
                    	 [replicate(cont_range[1], Nimg)] ] 
            print, "WARNING: Autorange with multiple contours !!"
        end

		log = arrayfy(log, Nimg)
		factor = arrayfy(factor, Nimg)
		colmap = arrayfy(colmap, Nimg)
		nlevels = arrayfy(nlevels, Nimg)
		cont_fact = arrayfy(cont_fact, Nimg)
		std_colmap = arrayfy(std_colmap, Nimg)
		cont_ext = arrayfy(cont_ext , Nimg)
		cont_log = arrayfy(cont_log , Nimg)
		cont_inv = arrayfy(cont_inv , Nimg)
		cont_thick = arrayfy(cont_thick , Nimg)
		cont_smooth = arrayfy(cont_smooth , Nimg)
		cont_factor = arrayfy(cont_factor , Nimg)
		cont_colmap = arrayfy(cont_colmap , Nimg)
		cont_charsize = arrayfy(cont_charsize , Nimg)
		cube = arrayfy(cube, Nimg)
		txtct = arrayfy(txtct, Nimg)
			
    end else begin ;Nimg = 1
        range = [[range[0]],[range[1]]]
        cont_range = [[cont_range[0]],[cont_range[1]]]
    end  

    ; init layout parameters
    doc_size = make_array(2,/double)
    vis_size = make_array(2,/double)
    img_size = make_array(2,/double)
    colbar_area_size = make_array(2,/double)
    colbar_size = make_array(2,/double)
    colbar_offset = make_array(2,/double)
    text_offset = make_array(2,/double)
    frame_offset =  make_array(2,/double)

    ; set layout
	if not keyword_set(XYimg) then $
	    XYimg = find_layout(Nimg)   
		; partition vis area in XYimg[0] x Xyimg[1] images
    
	if frame ne 0 then begin
        frame_offset[0] = 1.8                   ; [cm] offset for size annotations
        frame_offset[1] = 0.5
    end else begin
        frame_offset[0] = 0.0
        frame_offset[1] = 0.0
    end

    doc_size[0] = 21                            ; [cm], xsize of whole document
    vis_size[0] = doc_size[0] - frame_offset[0] ; [cm], size of visualisation area
    img_size[0] = vis_size[0]/XYimg[0]          ; [cm], size of single image in vis. area
    colbar_area_size[0] = vis_size[0]           ; [cm], size of area to reserve for colorbar

    img_size[1] = img_size[0]                   ; [cm], square images
    vis_size[1] = img_size[1]*XYimg[1]          ; [cm]
    colbar_area_size[1] = 0.25 * vis_size[1]     ; [cm]
    
    if Nimg gt 4 then $
        colbar_area_size[1] /= 2                ; Large images get a smaller colbar
    
    doc_size[1] = vis_size[1] + colbar_area_size[1] + frame_offset[1] ; [cm]
    
    colbar_size[0] = 0.8 * vis_size[0]                  
    colbar_size[1] = colbar_area_size[1]/4

    colbar_offset[0] = (vis_size[0]-colbar_size[0])/2   
    colbar_offset[1] = 0.5 * colbar_area_size[1]

    text_offset[0] = 0.05 * img_size[0]
    text_offset[1] = 0.1 * img_size[1]

    if keyword_set(nocolbar) then begin
        doc_size[1] -=  colbar_area_size[1]
        colbar_area_size[*] = 0
        colbar_size[*] = 0
        colbar_offset[*] = 0
    end
    charsize = doc_size[1] / 21           ;  norm to 21 cm page width
    ticklen = 0.03 * doc_size[1] / 21           ; norm to 21 cm page width

    if not keyword_set(frame_charsize) then $
        frame_charsize = min([0.8 , 0.6 * charsize * 3/XYimg[0]] )  

	; norm to 0.6 at 3 ximages

    ; set document
    set_plot, 'PS'
    !p.font = 2

    device, filename=fout,/color,/encapsulated,bits=24 $ 
        ,/portrait, xsize=doc_size[0],ysize=doc_size[1]

	loadct, 1, /silent
	polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=color(255) ; background white

    ; main loop
    img_pos = make_array(2,/double)

    for i=0, Nimg-1 do begin

        if not keyword_set(silent) then $
            print, 'Reading '+fin[i]

		if not keyword_set(cube) then begin

	        img = readfits(fin[i],ext=ext,/silent) * factor[i]

		end else begin 
			
		    img_cube = readfits(fin[i],ext=ext,/silent) * factor[i]

			img = img_cube[*,*,cube[i]]

		end

        if keyword_set(smooth_fwhm) then $  ; smoothing
            img = filter_image(img, fwhm=smooth_fwhm,/all_pixels)

        img = zoom_image(img, zoom)
            
        img = rebin(img, img_res, img_res)    ; change resolution
    
        this_range = range[i,*]

        img = set_range(img, range = this_range , silent=silent)
        
        colbar_range = this_range

        img = apply_scale(img, this_range, log[i], silent=silent)
        
        ; colormap
        if colmap[i] ge 0 then begin 

	        if std_colmap[i] eq 0 then $
    	        setcolor, colmap[i],silent=silent, /seq $
        	else $
            	loadct, colmap[i], silent=silent
		
		end else begin
            
			module_no =  head_find(fin[i], "Effect_Module", ext=ext)
            flag_no =  head_find(fin[i], "Effect_Flag",ext=ext)
            
            print, "Auto Color via Helix Table "+strn(module_no)+" "+strn(flag_no)

            helix_tables, module_no, flag_no, inv=invert_colmap, values=helix_values
		end

        ; output
        row = floor(i/float(XYimg[0]))
        column = (i - row*XYimg[0] )
        img_pos[0] = column * img_size[0] + frame_offset[0]
        img_pos[1] = (XYimg[1] - row - 1) * img_size[1] $
            + colbar_area_size[1] + frame_offset[1]
      
        tv, img, img_pos[0], img_pos[1], xsize=img_size[0] $
            , ysize=img_size[1],/centimeters
     
        ; annotation

 		set_txtct, txtct, i, img_res, img

        text_pos = make_array(2,/double)
        text_pos[0] = (img_pos[0] + text_offset[0]) / doc_size[0]  ; [normal]
        text_pos[1] = (img_pos[1] + img_size[1] - text_offset[1]) / doc_size[1]  ; [normal]
        
		if n_elements(text) gt 1 then begin

             xyouts, text_pos[0], text_pos[1], text[i], $
                col=color(txtct[i]), /normal, charsize=1.5*frame_charsize 
        
		end else if i eq 0 then begin 
        
			xyouts, text_pos[0], text_pos[1], text, $
                col=color(txtct[i]), /normal, charsize=1.5*frame_charsize
        end
        
        ; frame
        xstart = img_pos[0]/ doc_size[0]  ;[normal]
        xend =  (img_pos[0] + img_size[0]) / doc_size[0]
        dx = img_size[0] / doc_size[0]

        ystart = img_pos[1]/ doc_size[1]
        yend = (img_pos[1] + img_size[1]) / doc_size[1]
        dx = img_size[1] / doc_size[1]
    
        if i eq 0 then begin    ; find physical size

            xtitle='[kpc]'
            ytitle='[kpc]'

            xmax = head_find(fin[0],'XYSize', ext=ext)
            ymax = xmax
        end
    
        ytickname = replicate(' ',10) ; NULL ticknames
        xtickname = replicate(' ',10)

        if frame ne 0 then begin
            
			if column eq 0 then $
                plot, [1],/nodata,/noerase,position=[xstart,ystart,xend,yend]$
                    , xrange=[0,xmax], xtitle=' ',xstyle=1, xtickname=xtickname $
                    , yrange=[0,ymax], ytitle=ytitle, ystyle=1 , col=color(txtct[i])$
                    , charsize=frame_charsize,ticklen=ticklen, xthick=xythick, ythick=xythick
        
            if  row eq XYimg[1]-1 then $
                 plot, [1],/nodata,/noerase,position=[xstart,ystart,xend,yend]$
                    , xrange=[0,xmax], xtitle=xtitle,xstyle=1 , col=color(txtct[i])$
                    , yrange=[0,ymax], ytitle=' ', ystyle=1, ytickname=ytickname $
                    , charsize=frame_charsize,ticklen=ticklen, xthick=xythick, ythick=xythick

            plot, [1],/nodata,/noerase,position=[xstart,ystart,xend,yend] , col=color(txtct[i])$
                , xrange=[0,xmax], xtitle=' ',xstyle=1, xtickname=xtickname $
                , yrange=[0,ymax], ytitle=' ', ystyle=1, ytickname=ytickname $
                , charsize=frame_charsize,ticklen=ticklen, xthick=xythick, ythick=xythick

        end else begin      ; annotation

 			set_txtct, txtct, i, img_res, img
        
            text_pos = make_array(2,/double)
            text_pos[0] = (img_pos[0] + 0.7 * img_size[0]  ) / doc_size[0]  ; [normal]
            text_pos[1] = (img_pos[1] + 0.1 * img_size[1]) / doc_size[1]  ; [normal]

            if i eq 0 then begin

                xysize = head_find(fin[0],'XYSize', ext=ext) / double(zoom)
    
                if not exist(bartext) then begin
                    bartext = '['+strn(xysize, len=3)+' kpc]!U2!N'
                    
				if xysize gt 1e3 then $
                    bartext = '['+strn(xysize/1000L,len=3)+' Mpc]!U2!N'
                end

                xyouts, text_pos[0], text_pos[1], bartext, $
                    col=color(txtct[i]), /normal, charsize=1.5*frame_charsize
            end

            plot, [1],/nodata,/noerase,position=[xstart,ystart,xend,yend]$
                , xrange=[0,xmax], xtitle=' ',xstyle=1, xtickname=replicate(' ',10) $
                , yrange=[0,ymax], ytitle=' ', ystyle=1, ytickname=replicate(' ',10)  $
                , charsize=0,ticklen=ticklen, xthick=xythick, ythick=xythick $
                , col=color(txtct[i]), xminor=frame_nminor, yminor=frame_nminor
        end

		; contour plot
    	if keyword_set(cont_fin) then begin ; needs coyote lib

	        if not keyword_set(cont_thick) then $
        	    cont_thick = 4

			if not keyword_set(cont_label) then $
				cont_label = replicate('', Nimg)

    		if not keyword_set(cont_charsize[i]) then $
				cont_charsize[i] = charsize / 5

 			if cont_inv[i] eq 0 then $
				cont_col = replicate('White', Nimg) $
			else $
				cont_col = replicate('Black', Nimg) 
			
			print, 'Reading Contours from ', cont_fin[i]

			img = readfits(cont_fin[i], /silent, ext=ext) * cont_factor[i]

	        img = zoom_image(img, zoom)

        	img = rebin(img, img_res, img_res)    ; change resolution

	        if  keyword_set(cont_smooth) then $
    	      img = gauss_smooth(img, cont_smooth[i] )
        
			this_range = cont_range[i,*]
    
    	    img = set_range(img, range=this_range, silent=silent)

	        if not keyword_set(cont_levels) then $
        	    cont_levels = make_contour_levels(this_range, nlevels[i], $
					cont_log[i]) $
    	    else $
	            nLevels = n_elements(cont_levels)

        	cgcontour, img, /noerase,col=cont_colmap[i],levels=cont_levels, $ 
    	        position=[xstart,ystart,xend,yend], $
				xtickname=replicate(' ',10), $
	            c_thick=replicate(cont_thick[i], nlevels), $
				xthick=xythick, ythick=xythick, $
            	c_charsize=cont_charsize[i], ytickname=replicate(' ',10), $
        	    xticks=1, yticks=1, xminor=1, yminor=1, $ 
				c_colors=cont_col[i], label=cont_label[i] 
    	end

    end

    ; colorbar
    if not keyword_set(nocolbar) then begin

        ystart = colbar_offset[1]
        yend = colbar_offset[1] + colbar_size[1]
        dy = colbar_size[1]

        xstart = frame_offset[0] + colbar_offset[0]
        xend = frame_offset[0] + colbar_offset[0] + colbar_size[0]
        dx = colbar_size[0]    

        pos = [xstart,ystart,xend,yend]

        bar = BINDGEN(256) # REPLICATE(1B, barthick)

        tvscl, bar, xstart, ystart , xsize=dx, ysize=dy,/centimeters

        loadct, 0,/silent

        if not keyword_set(unit) then begin ; try get unit from header of img[0]
            str1 = head_find(fin[0], 'Module Name',ext=ext)
            str2 = head_find(fin[0], 'Unit',ext=ext)

            if (size(str1))[1] ne 7 or (size(str2))[1] ne 7 then $ ;info not found
                unit='[arb. units]' $
            else $
                unit = strtrim(str1,2)+' '+strtrim(str2,2)


        end
        
        if  not keyword_set(silent) then $
            print, 'Found unit '+unit

        xstart /= doc_size[0] ;convert to [normal]
        xend /= doc_size[0]

        ystart /= doc_size[1]
        yend /= doc_size[1]

        plot, [1.1],/nodata,/noerase,position=[xstart,ystart,xend,yend], $
            xrange=colbar_range, xlog=log[0],yticks=1, xtitle=unit,xstyle=1, $
            yminor=1,xticklen=0.28, ytickname=REPLICATE(' ',2), $
            charsize = 1.5*frame_charsize, xthick=xythick, ythick=xythick, xticks=bar_nticks
    end



; overplot polarisation vectors
    if keyword_set(vect) then begin
        res = 64

        ang = readfits(fin,ext=4,/silent)
        Ipol = readfits(fin, ext=3,/silent)

        u = cos(ang+!pi/2)
        v = sin(ang+!pi/2)
      
        u = rebin(u, res, res)
        v = rebin(v, res, res)
        Ipol = rebin(Ipol, res, res)

        ipol = filter_image(Ipol, fwhm=res/10)

        u *= Ipol
        v *= Ipol
        
        l = sqrt(u*u + v*v)

        bad = where(l LT 0.01 * max(l),cnt)
        
        if cnt NE 0 then begin
            u[bad] = 1d35
            v[bad] = 1d35
        end

        vecplot, u,v,xticks=1,yticks=1,ytickname=REPLICATE(' ',2) $
            ,xtickname=REPLICATE(' ',2),xticklen=0.01,/noerase $
            ,position=[0,0.2,1,1],thick=0.5, length=1
    end
    
    device, /close

    if keyword_set(pdf) then $
		spawn, 'epstopdf '+fout+' && rm '+fout

    return
end


