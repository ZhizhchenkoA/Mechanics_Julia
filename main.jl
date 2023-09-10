### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ f654666c-36bb-4854-ae63-3f0a1bb0343f
begin
	using Plots
	using ProgressBars
end

# ╔═╡ d3f35948-a39c-462c-9412-8bbefa334fe2
html"""
	<h1 align="center">Научно-исследовательская работа</h1>
	
	<h2 align="center"><b>"Численное моделирование механических систем"</b></h2>
	"""

# ╔═╡ a3e811d7-0a06-44a6-a9db-0d21983b212a
begin
	md"""
	_Выполнили:_
	* Липская Мария
	* Карасева Ольга
	* Жижченко Александр

	
	Научный руководитель:
	* Федоров Владимир Сергеевич
	
	*Цели работы*:
	* научиться численно моделировать механические системы
	* изучить решения ограниченной задачи трёх тел
	* разобраться с положением равновесия в точках Лагранжа
	"""
end

# ╔═╡ 56257e88-88f1-4c52-a9b4-382f6de339f6
begin
	mu = 0.2
	function u(x, y)
	    """Потенциал тела U(x, y)"""
	    r1 = sqrt((x + mu)^2 + y^2)
	    r2 = sqrt((x - 1 + mu)^2 + y^2)
	    return (x^2 + y^2) / 2 + (1 - mu) / r1 + mu / r2 + ((1 - mu) * mu) / 2
	end
	
	function h(x, y, px, py)
	    """Гамильтониан системы"""
	    return ((px + y)^2 + (py - x)^2) / 2 - u(x, y)
	end
	
	
	function diff_u_x(x, y)
	    """Производная в точке от U(x, y) по x"""
	    delta = 1e-5
	    return  - (u(x + delta / 2, y) - u(x - delta / 2, y)) / delta
	end
	
	function diff_u_y(x, y)
	    """Производная в точке от U(x, y) по y"""
	    delta = 1e-5
	    return  - (u(x, y + delta / 2) - u(x, y - delta / 2)) / delta
	end
	
	function x_dot(x, y, px, py)
	    """Скорость по x"""
	    return px + y
	end
	
	function y_dot(x, y, px, py)
	    """Скорость по y"""
	    return py - x
	end
	
	function px_dot(x, y, px, py)
	    """Производная по px"""
	    return py - x + diff_u_x(x, y)
	end
	
	function py_dot(x, y, px, py)
	    """Производная по py"""
	    return -px - y + diff_u_y(x, y)
	end
	
	function x_ddot(x, y, px, py)
	    """Вторая производная по x """
	    return px_dot(x, y, px, py) + y_dot(x, y, px, py)
	end
	
	function y_ddot(x, y, px, py)
	    """Вторая производная по Y """
	    return py_dot(x, y, px, py) - x_dot(x, y, px, py)
	end
	
	function all_dots(x, y, px, py, with_second=true)
	    """Все производные: """
	    if with_second
	        return (x_dot(x, y, px, py), y_dot(x, y, px, py), px_dot(x, y, px, py), py_dot(x, y, px, py), x_ddot(x, y, px, py), y_ddot(x, y, px, py))
	    else
	        return (x_dot(x, y, px, py), y_dot(x, y, px, py), px_dot(x, y, px, py), py_dot(x, y, px, py))
	    end
	end
	
	function velocities_half_delta(x, y, px, py, delta_t)
	    x_new = x - x_dot(x, y, px, py) * delta_t / 2
	    y_new = y - y_dot(x, y, px, py) * delta_t / 2
	    px_new = px - px_dot(x, y, px, py) * delta_t / 2
	    py_new = py - py_dot(x, y, px, py) * delta_t / 2
	    return (x_dot(x_new, y_new, px_new, py_new), y_dot(x_new, y_new, px_new, py_new))
	end
	
	function W(x::Float64, y::Float64)
	    """Функция нормализации значения"""
	    return -log((abs(diff_u_x(x, y)) + abs(diff_u_y(x, y))), ℯ)
	end
	
	function h_const(x, y, px=0, py=0)
	    """Подсчёт постоянной Якоби"""
	    return 0.5 * (x_dot(x, y, px, py) ^ 2 + y_dot(x, y, px, py) ^ 2) - u(x, y)
	end
	
	function puancare(X, Y, delta_t)
		
	    """Построение сечения Пуанкаре по заданной траектории"""
		x_result = []
		x_dot_result = []
	    for (n, y) in enumerate(Y[1:length(Y) - 1])
	        
	        if y * Y[n + 1] <= 0
	            
	            append!(x_result, X[n])
	            append!(x_dot_result, (X[n + 1] - X[n]) / delta_t)
	        end
	    end
	    return (x_result, x_dot_result)
	end
	
end

# ╔═╡ 8ff68780-6554-4873-9abb-aa118e116f62


# ╔═╡ 7befb0a2-e81f-4d8b-bd09-4f253dc36325
begin
	lim_x_min, lim_x_max = -2, 2
	step_x = 0.01
	lim_y_min, lim_y_max = -2, 2
	step_y = 0.01
	DAT = zeros(length(lim_x_min:step_x:lim_x_max), length(lim_y_min:step_y:lim_y_max))
	
	for (j, y) in enumerate(lim_x_min:step_x:lim_x_max)
	    for (i, x) in enumerate(lim_y_min:step_y:lim_y_max)
	        DAT[j, i] = W(x, y)
	    end
	end
	
	
end

# ╔═╡ 1e8db92f-d780-4052-9381-87d0719659cb
html"""
<h1 align="center"> Эпизод II<br> <b>Нахождение координат точек Лагранжа</b></h1>
<hr>
"""

# ╔═╡ ad124e1c-e3bf-47f7-b7a5-f1d97e00a14c
md"""

Для нахождения координат точек Лагранжа, нами использовались массивы,
```python
mx = np.arange(lim_x_min, lim_x_max, step_x)
my = np.arange(lim_y_min, lim_y_max, step_y)
```
которые задаются минимальными значениями, максимальными значениями и шагом изменения с помощью функции `numpy.arrange`. В итоговый массивы точек Лагранжа `ans_x` и `ans_y` заносились координаты, в которых $\frac{∂}{\partial x}U(x, y) = 0$ и $\frac{∂}{\partial y}U(x, y) = 0$ (реализовано с помощью приращений и теоремы Лагранжа)
Из массива итоговых точек были исключены координаты центра масс и тел, и получившиеся координаты были добавлены в массивы `L_x` и `L_y` и отсортированы по номеру точки.  
Координаты этих точек также задаются формулами:

$\begin{gather}
L_1 \:(R(1 - \sqrt[3]{\frac{\mu}{3}},\; 0))\\
L_2 \:(R(1 + \sqrt[3]{\frac{\mu}{3}}, \; 0)\\
L_3 \:(R(1 - \frac{5}{12}\mu, \; 0)\\
L_4 \:(\frac{R}{2} \cdot \frac{m_1-m_2}{m_1+m_2}, \; \frac{\sqrt[3]{3}}{2}R)\\
L_5 \:(\frac{R}{2} \cdot \frac{m_1-m_2}{m_1+m_2}, \; -\frac{\sqrt[3]{3}}{2}R)\\
\end{gather}$

где:  

- ㅤ$m_1$ масса большего тела  

- ㅤ$m_2$ масса меньшего тела

- ㅤ$R$ расстояние между телами
"""

# ╔═╡ e41e33a5-212f-403f-b51c-ae94dfab293b
begin
	ans_x = []
	ans_y = []
	epsilon = 0
	
	for x in lim_x_min:step_x:lim_x_max
	    for y in lim_y_min:step_y:lim_y_max
	         if (u(x, y) - u(x - step_x, y)) * (u(x + step_x, y) - u(x, y)) < epsilon && (u(x, y) - u(x, y - step_y)) * (u(x, y + step_y) - u(x, y)) < epsilon
	                append!(ans_x, x)
	                append!(ans_y, y)
	        end
	    end
	end
	L_x = [ans_x[7], ans_x[9], ans_x[1], ans_x[4], ans_x[3]]
	L_y = [ans_y[7], ans_y[9], ans_y[1], ans_y[4], ans_y[3]]
	println(L_x)
	println(L_y)
end

# ╔═╡ cfec2969-c486-43e8-aeed-f53157a55563
html"""
<h1 align="center"> Эпизод III<br> <b>Моделирование движения вблизи точек Лагранжа</b></h1><hr>
"""

# ╔═╡ 9a3aedf3-4958-4284-81f5-654f7835fae2
md"""

Движение точки в задаче трёх тел описывается уравнениями движения, описанными во введении. В программе перемещение точки задаётся функцией `move(x, y, px=0, py=0, t=100000, scfale_large)`, где `x` и `y` - обязательные параметры, являющиеся начальными координатами, а `px` и `py` начальные импульсы тела, `t` - общее время движения тела, `scale_large` (если `True`, то показано движение тела относительно системы тел, если `False`, то демонстрируется движение у точки Лагранжа)
"""

# ╔═╡ de3a78b5-bcd0-4537-8354-1418f449be4a
function move_leapfrog(x::Float64, y::Float64; px::Float64=0.0, py::Float64=0.0, t::Int64=100000)
    delta_t = 1e-4
    x_vel, y_vel = velocities_half_delta(x, y, px, py, delta_t)
    x_arr = [x]
    y_arr = [y]
    for i in 0:1:t
        derivatives = all_dots(x_arr[length(x_arr)], y_arr[length(y_arr)], px, py)

        x_vel += derivatives[5] * delta_t
        y_vel += derivatives[6] * delta_t

        push!(x_arr, x_arr[length(x_arr)] + x_vel * delta_t)
        push!(y_arr, y_arr[length(y_arr)] + y_vel * delta_t)

        px += derivatives[3] * delta_t
        py += derivatives[4] * delta_t
    end
    return (x_arr, y_arr)
end

# ╔═╡ 911f528b-816a-4c6e-9538-2530309f70c6
function move_leapfrog_gif(x::Float64, y::Float64; px::Float64=0.0, py::Float64=0.0, t::Int64=100000, frames::Int64=100)
    delta_t = 1e-4
    x_vel, y_vel = velocities_half_delta(x, y, px, py, delta_t)
    x_arr = [x]
    y_arr = [y]
	#scatter([-mu, 1-mu], [0, 0])
    @gif for i in tqdm(0:1:t)
        derivatives = all_dots(x_arr[length(x_arr)], y_arr[length(y_arr)], px, py)

        x_vel += derivatives[5] * delta_t
        y_vel += derivatives[6] * delta_t

        push!(x_arr, x_arr[length(x_arr)] + x_vel * delta_t)
        push!(y_arr, y_arr[length(y_arr)] + y_vel * delta_t)

        px += derivatives[3] * delta_t
        py += derivatives[4] * delta_t
		plot(x_arr, y_arr)
	end every frames
    
end

# ╔═╡ f06cf227-88de-4f3f-9fc5-1714330a8e3a
begin
	L1 = move_leapfrog(L_x[1] - 1e-6, 0.0, px=0.0, py=0.4371024, t=14520)
	Plots.plot(L1[1], L1[2])
end

# ╔═╡ 32629dd1-a7fd-45d9-9623-d5665586f489
begin
	L2 = move_leapfrog(L_x[2], L_y[2], px=0.0, py=1.27183, t=21500)
	Plots.plot(L2[1], L2[2])
end

# ╔═╡ b8411d90-48e6-403c-ad68-d6671bb53dd6
begin
	L3 = move_leapfrog(L_x[3], L_y[3], px=0.0, py=-1.0855, t=23500)
	Plots.plot(L3[1], L3[2])
end

# ╔═╡ 0f50bb6b-55ca-446a-9ee3-18b528b947d4
begin
	L4 = move_leapfrog(L_x[4], L_y[4], px=-1.0, py=0.0, t=1000000)
	Plots.plot(L4[1], L4[2])
end

# ╔═╡ ec834d95-fba6-4b7f-9ab4-804113dabf1b
begin
	L5 = move_leapfrog(L_x[5], L_y[5], px=0.9, py=0.5, t=1000000)
	Plots.plot(L5[1], L5[2])
end

# ╔═╡ 93267c29-cc40-4515-a1fb-28c806d4cbb5
begin
	X, Y = move_leapfrog(L_x[1], L_y[1], px=0.0, py=0.0, t=300000)
	Plots.plot(X, Y)
	Plots.scatter!([-mu, 1-mu], [0, 0])
end

# ╔═╡ bb78e7c5-53bf-4c3d-8fec-a9f5a6adea5a
begin
	delta_t = 1e-4
	
	ans_x1, ans_y1 = puancare(X, Y, delta_t)
	Plots.scatter(ans_x1, ans_y1)
end

# ╔═╡ bd0514de-cc8d-485a-b0e5-46c8d130988f
begin
	x = 0.4  # @param {type:"slider", min:-1, max:1, step:0.001}
	y = 0.4  # @param {type:"slider", min:-1, max:1, step:0.001}
	px = 0.334  # @param {type:"slider", min:-1, max:1, step:0.001}
	py = -0.5  # @param {type:"slider", min:-1, max:1, step:0.001}
	t = 550000
	X1, Y1 = move_leapfrog(x, y, px=px, py=py, t=t)
	
	Plots.plot(X1, Y1)
	Plots.scatter!([-mu, 1-mu], [0, 0])
end

# ╔═╡ d23baf85-ce21-4a93-bc95-a6643cd59ce7
begin
	ans_x2, ans_y2 = puancare(X1, Y1, 1e-4)
	Plots.scatter(ans_x2, ans_y2)
end

# ╔═╡ 8eb16823-767b-4fcf-98da-3da540ca66fd
md"""
__Вывод интеграла Якоби, единственной формы энергии, сохраняющейся в ограниченной задаче трёх тел:__

Интегрируем выражение 

$\frac{dU(x, y)}{dt} = \dot xU_x + \dot yU_y = \dot x\ddot x + \dot y\ddot y$

и получаем

$\begin{gather}
2U = \dot x^2 + \dot y^2 + C\\
C = -(\dot x^2 + \dot y^2) + 2U
\end{gather}$

Поделим это выражение на 2

$C_j = \frac{(\dot x^2 + \dot y^2)}{2} - U$

где $C_j$ - постоянная Якоби.

Эквипотенциальные поверхности задаются с помощью следующего отображения

$H(h) = \{ (x, y) \:\colon U(x, y) + C_j \}$

В нашем случае нам необходимо оценить границы, в пределах которых движется третье тело, поэтому за начальную точку постоянной Якоби $C_j$ была выбрана точка Лагранжа $L_1$, импульсы $p_x$ и $p_y$.
Необходимые нам граничные точки удовлетворяют условию: $C_j + U(x, y) = 0$


 Код ниже вычисляет постоянную Якоби и на основании него строит эквипотенциальную поверхность (функция `h_const(x, y, px, py)`). С отклонением `epsilon` значение функции $U(x, y) + C_j = 0$
 """

# ╔═╡ 1207920a-0ebc-46de-b426-3abea1664eef
begin
	ans_xx, ans_yy = [], []
	ϵ = 0.01
	const_h = h_const(L_x[1], L_y[1], 0, 0)
	for i in -2.0:0.01:2.0
	    for j in -2.0:0.01:2.0
	        if abs(u(i, j) + const_h) <= ϵ
	            append!(ans_xx, i)
	            append!(ans_yy, j)
	        end
	    end
	end
	a, b = move_leapfrog(L_x[1], 0.4, px=0.0, py=0.0, t=10000000)
	scatter(ans_xx, ans_yy)
	plot!(a, b)
end

# ╔═╡ 30358f2c-0329-4b59-82e4-f208e8596f25
md"""
Тест при большом расстоянии третьего тела от объектов 
"""

# ╔═╡ c361a162-30ee-47b8-88f4-2e59225da70c
begin
	ɼ1, ɼ2 = move_leapfrog(1000.0, 10000000.0, px=1000000.0, py=10000000.0, t=1000000)
	plot(ɼ1, ɼ2)
end

# ╔═╡ 56d02fea-a060-4634-abe9-082c706fc0f3
md"""
Построение анимации в сидерической системе отсчёта
"""

# ╔═╡ 1c0ad33c-2490-4922-9c61-8b5595cef356
@fastmath function animate_non_rotating(x::Float64, y::Float64; px::Float64=0.0, py::Float64=0.0, t::Int64=100000, frames::Int64 = 200)
	ro1, ro2 = -mu, 1 - mu
	psi = 0
	circle_dots = range(0, 2π, length = 100)

		
	delta_t = 1e-4
    x_vel, y_vel = velocities_half_delta(x, y, px, py, delta_t)
    x_arr = [x]
    y_arr = [y]
	points = [-mu, 1-mu, x], [0, 0, y]
	scatter(points)
	
	anim = @animate for i in 0:1:t
        derivatives = all_dots(x_arr[length(x_arr)], y_arr[length(y_arr)], px, py)

        x_vel += derivatives[5] * delta_t
        y_vel += derivatives[6] * delta_t

        push!(x_arr, x_arr[length(x_arr)] + x_vel * delta_t)
        push!(y_arr, y_arr[length(y_arr)] + y_vel * delta_t)
		x, y = x_arr[length(x_arr)], y_arr[length(y_arr)]
        px += derivatives[3] * delta_t
        py += derivatives[4] * delta_t

		psi += delta_t
		x1, y1 = ro1 * cos(psi), ro1 * sin(psi)
		x2, y2 = ro2 * cos(psi), ro2 * sin(psi)
		if i % frames == 0

			
			x_1 = x * cos(psi) - y * sin(psi)
			y_1 = x * sin(psi) + y * cos(psi)
			scatter([x1, x2], [y1, y2], xlim=(-2, 2), ylim=(-2, 2), color="black")
			scatter!([x_1], [y_1], color="magenta")
			plot!(@. sin(circle_dots) * ro1, @. cos(circle_dots) * ro1) 
			plot!(@. sin(circle_dots) * ro2, @. cos(circle_dots) * ro2) 
		end
	end every frames
	gif(anim, "main.gif", fps=60)
end

# ╔═╡ 5939245c-0571-4e24-a39b-53d813a61e62
animate_non_rotating(L_x[5], L_y[5], px=0.0, py=-1.0, t=100000, frames=200)


# ╔═╡ Cell order:
# ╠═f654666c-36bb-4854-ae63-3f0a1bb0343f
# ╟─d3f35948-a39c-462c-9412-8bbefa334fe2
# ╟─a3e811d7-0a06-44a6-a9db-0d21983b212a
# ╠═56257e88-88f1-4c52-a9b4-382f6de339f6
# ╠═8ff68780-6554-4873-9abb-aa118e116f62
# ╠═7befb0a2-e81f-4d8b-bd09-4f253dc36325
# ╟─1e8db92f-d780-4052-9381-87d0719659cb
# ╟─ad124e1c-e3bf-47f7-b7a5-f1d97e00a14c
# ╠═e41e33a5-212f-403f-b51c-ae94dfab293b
# ╟─cfec2969-c486-43e8-aeed-f53157a55563
# ╟─9a3aedf3-4958-4284-81f5-654f7835fae2
# ╠═de3a78b5-bcd0-4537-8354-1418f449be4a
# ╠═911f528b-816a-4c6e-9538-2530309f70c6
# ╠═f06cf227-88de-4f3f-9fc5-1714330a8e3a
# ╠═32629dd1-a7fd-45d9-9623-d5665586f489
# ╠═b8411d90-48e6-403c-ad68-d6671bb53dd6
# ╠═0f50bb6b-55ca-446a-9ee3-18b528b947d4
# ╠═ec834d95-fba6-4b7f-9ab4-804113dabf1b
# ╠═93267c29-cc40-4515-a1fb-28c806d4cbb5
# ╠═bb78e7c5-53bf-4c3d-8fec-a9f5a6adea5a
# ╠═bd0514de-cc8d-485a-b0e5-46c8d130988f
# ╠═d23baf85-ce21-4a93-bc95-a6643cd59ce7
# ╟─8eb16823-767b-4fcf-98da-3da540ca66fd
# ╠═1207920a-0ebc-46de-b426-3abea1664eef
# ╟─30358f2c-0329-4b59-82e4-f208e8596f25
# ╠═c361a162-30ee-47b8-88f4-2e59225da70c
# ╟─56d02fea-a060-4634-abe9-082c706fc0f3
# ╠═1c0ad33c-2490-4922-9c61-8b5595cef356
# ╠═5939245c-0571-4e24-a39b-53d813a61e62
