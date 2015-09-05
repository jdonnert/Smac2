	pro helix_tables, module, flag, inv=inv, values=values

    nModule = 100
    nFlag = 100
    
    ; [ start, rots, hue, gamma]
    setup = make_array(4, nModule, nFlag, val=0, /float)

    ;std values
    for i = 0,nModule-1 do $
	for j = 0,nFlag-1 do $
    	   setup[*,i,j] = [ -0, 0, 2.5, 1 ]

    ; density
    setup[*, 0, 0] = [ 3, 3, 3, 2 ]

    ; velocities
    setup[*, 1, 0] = [ 1, 1.5, 2, 1 ]
    
    ; x-rays
    setup[*, 2, 0] = [ 2, -2, 3, 1 ]

    ; Temperature
    setup[*, 4, 0] = [ 0.9, 0, 1, 1 ] ; Mass Weighted
    setup[*, 4, 1] = [ 1, 0, 1, 1 ] ; Sound Speed
    setup[*, 4, 2] = [ 1, 0, 1, 1 ] ; Emission Weighted
    setup[*, 4, 3] = [ 1.3, -0.5, 2, 2 ] ; Spectroscopic

    ; Pressure
    setup[*, 5, 0] = [ 0, 3, 1, 1 ]

    ; magnetic field
    setup[*, 6, 0] = [ 3, 0, 3, 2]

    ; Compton-Y / SZ
    setup[*, 7, 0] = [ 3, -1, 4, 1 ]
    setup[*, 7, 1] = [ 2, 1, 4, 1 ]

    ; Dm density
    setup[*, 10, 0] = [ 3, -1, 4, 1 ]
    setup[*, 11, 0] = [ 1, 0, 1.4, 1 ]

    if keyword_set(values) then $
	setup[*,module,flag] = values

    ; set
    setchelix, setup[0,module, flag], setup[1, module, flag], $
        setup[2,module, flag], setup[3, module, flag], inv=inv

    return
end
