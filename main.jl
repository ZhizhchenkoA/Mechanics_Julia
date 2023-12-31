### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ f654666c-36bb-4854-ae63-3f0a1bb0343f
begin
	using Plots
	using ProgressBars
	using Parameters
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
	\
	"""
end

# ╔═╡ f997a055-a870-430a-8aca-5060c2fd16cb
begin
	@with_kw mutable struct Body
		x::Array{Float64}= []
		y::Array{Float64} = []
		px::Float64 = 0
		
		py::Float64 = 0
		mu::Float64 = 0.2

	
		function Body(x, y, px=0, py=0, mu=0.2)
			new([x], [y], px, py, mu)
		end
	end

	function change_body!(body::Body, x::Float64, y::Float64, px::Float64, py::Float64)
		push!(body.x, x)
		push!(body.y, y)
		body.px = px
		body.py = py
		nothing
	end

	function body_cords(body::Body)
		return [body.x[end], body.y[end], body.px, body.py, body.mu]
	end
end

# ╔═╡ d3f8fd38-ba29-4234-a316-724cfb5ff183
struct Point
	x::Float64
	y::Float64
end

# ╔═╡ 56257e88-88f1-4c52-a9b4-382f6de339f6
begin
	function u(x::Float64, y::Float64, mu::Float64=0.2)
	    """Потенциал тела U(x, y)"""
	    r1 = sqrt((x + mu)^2 + y^2)
	    r2 = sqrt((x - 1 + mu)^2 + y^2)
	    return - ((x^2 + y^2) / 2 + (1 - mu) / r1 + mu / r2 + ((1 - mu) * mu) / 2)
	end
	
	function h(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Гамильтониан системы"""
	    return ((px + y)^2 + (py - x)^2) / 2 - u(x, y, mu)
	end
	
	
	function diff_u_x(x::Float64, y::Float64, mu::Float64=0.2)
	    """Производная в точке от U(x, y) по x"""
	    delta = 1e-10
	    return  (u(x + delta / 2, y, mu) - u(x - delta / 2, y, mu)) / delta
	end
	
	function diff_u_y(x::Float64, y::Float64, mu::Float64=0.2)
	    """Производная в точке от U(x, y) по y"""
	    delta = 1e-10
	    return  (u(x, y + delta / 2) - u(x, y - delta / 2)) / delta
	end


	function x_dot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Скорость по x"""
	    return px + y
	end
	
	function y_dot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Скорость по y"""
	    return py - x
	end
	
	function px_dot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Производная по px"""
	    return py - x + diff_u_x(x, y, mu)
	end
	
	function py_dot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Производная по py"""
	    return -px - y + diff_u_y(x, y, mu)
	end
	
	function x_ddot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Вторая производная по x """
	    return px_dot(x, y, px, py, mu) + y_dot(x, y, px, py, mu)
	end
	
	function y_ddot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
	    """Вторая производная по Y """
	    return py_dot(x, y, px, py) - x_dot(x, y, px, py)
	end

	function px_ddot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
		
		return py_dot(x, y, px, py, mu) - x_dot(x, y, px, py, mu) - diff_u_ddx(x, y, mu)

	end

	function py_ddot(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
		return -px_dot(x, y, px, py, mu) - y_dot(x, y, px, py, mu) - diff_u_ddy(x, y, mu)

	end
	
	function all_dots(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2; with_second=true, with_impulses=false)
	    """Все производные: """
	    if with_second && with_impulses
	        return (x_dot(x, y, px, py, mu), y_dot(x, y, px, py, mu), px_dot(x, y, px, py, mu), py_dot(x, y, px, py, mu), x_ddot(x, y, px, py, mu), y_ddot(x, y, px, py, mu), px_ddot(x, y, px, py, mu), py_ddot(x, y, px, py, mu))
		elseif with_second
			return (x_dot(x, y, px, py, mu), y_dot(x, y, px, py, mu), px_dot(x, y, px, py, mu), py_dot(x, y, px, py, mu), x_ddot(x, y, px, py, mu), y_ddot(x, y, px, py, mu))
	    else
	        return (x_dot(x, y, px, py, mu), y_dot(x, y, px, py, mu), px_dot(x, y, px, py, mu), py_dot(x, y, px, py, mu))
	    end
	end
	
	function q_dots(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
		return [x_dot(x, y, px, py, mu), y_dot(x, y, px, py, mu), 0, 0]
	end

	function p_dots(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2)
		return [0, 0, px_dot(x, y, px, py, mu), py_dot(x, y, px, py, mu)]
	end

	
	function initial_velocities(x::Float64, y::Float64, px::Float64, py::Float64, mu::Float64=0.2, delta_t::Float64=1e-4)
	    x_new = x + x_dot(x, y, px, py, mu) * delta_t
	    y_new = y + y_dot(x, y, px, py, mu) * delta_t 
	    px_new = px + px_dot(x, y, px, py, mu) * delta_t
	    py_new = py + py_dot(x, y, px, py, mu) * delta_t 
	    return (x_dot(x_new, y_new, px_new, py_new, mu), y_dot(x_new, y_new, px_new, py_new, mu), px_dot(x_new, y_new, px_new, py_new, mu), py_dot(x_new, y_new, px_new, py_new, mu))
	end

	function W(x::Float64, y::Float64, mu::Float64=0.2)
	    """Функция нормализации значения"""
	    return -log((abs(diff_u_x(x, y, mu)) + abs(diff_u_y(x, y, mu))), 10)
	end
	
	function h_const(x::Float64, y::Float64, px::Float64=0, py::Float64=0, mu::Float64=0.2)
	    """Подсчёт постоянной Якоби"""
	    return - (x_dot(x, y, px, py, mu) ^ 2 + y_dot(x, y, px, py, mu) ^ 2) + 2 * u(x, y, mu)
	end
	
	function puancare(X, Y, delta_t)
		
	    """Построение сечения Пуанкаре по заданной траектории"""
		x_result = []
		x_dot_result = []
	    for (n, y) in enumerate(Y[1:end-1])
	        
	        if y * Y[n + 1] <= 0
	            
	            append!(x_result, X[n])
	            append!(x_dot_result, (X[n + 1] - X[n]) / delta_t)
	        end
	    end
	    return (x_result, x_dot_result)
	end
	nothing
end

# ╔═╡ c1590f21-e4eb-48cb-8c0d-fb98fbfed059
function gradient(x::Float64, mu::Float64=0.2; iterations::Int64=100)
	λ = 0.01
	for i ∈ 1:1:100
		x += λ*diff_u_x(x, 0.0, mu)
	end
	return x
end

# ╔═╡ fc385955-4924-4326-ab8e-d37ecd5f65e4
@with_kw mutable struct Lagrange
	L1::Point
	L2::Point
	L3::Point
	L4::Point
	L5::Point

	function Lagrange(mu::Float64)
		l1::Point = Point(gradient(1 - ((mu/ (1 - mu)) / 3)^(1/3), mu), 0)
		l2::Point = Point(gradient(1 + ((mu/ (1 - mu)) / 3)^(1/3), mu), 0)
		l3::Point = Point(gradient(- 1 - 5*(mu/ (1 - mu)) / 12, mu), 0)
		l4::Point = Point(1/2 - mu, sqrt(3) / 2)
		l5::Point = Point(1/2 - mu, sqrt(3) / 2)
		new(l1, l2, l3, l4, l5)
	end
end

# ╔═╡ 41788216-b144-4460-bded-e87992c3d57a
Lagrange(0.2).L4

# ╔═╡ 9f7c73f7-226e-4c04-9741-834f46531deb
html"""
<h1 align="center"> Эпизод I<br> <b>Визуализация точек Лагранжа, полости Роша и сферы Хилла</b></h1>
"""

# ╔═╡ 8ff68780-6554-4873-9abb-aa118e116f62
md"""
Для наглядного представления точек Лагранжа нами была выбрана функция `imshow` из пакета `matplotlib.pyplot`. Она задаёт изображение, как двумерный массив пикселей (в нашем коде массив `DAT`). Его заполнение осуществлялось с помощью перебора значений `i` от 0 до `(lim_x_max - lim_x_min) / step_x` и `j` от 0 до `(lim_y_max - lim_y_min) / step_y`. Каждому пикселю `DAT[i][j]` присваивалось значение `W(i, j)`. Для большей наглядности результата изображение было выведено в чёрно-белом формате `cmap='gray'`
"""

# ╔═╡ 2a4aa267-796c-4090-ad20-3ce806868726
begin
	lim_x_min, lim_x_max = -2, 2
		step_x = 0.01
		lim_y_min, lim_y_max = -2, 2
		step_y = 0.01
end

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

# ╔═╡ 7befb0a2-e81f-4d8b-bd09-4f253dc36325
begin
	
	DAT = zeros(length(lim_x_min:step_x:lim_x_max), length(lim_y_min:step_y:lim_y_max))

	for (j, y) in enumerate(lim_x_min:step_x:lim_x_max)
	    for (i, x) in enumerate(lim_y_min:step_y:lim_y_max)
	        DAT[j, i] = W(x, y)
	    end
	end
	print(DAT[1])
	heatmap(DAT)
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
L_1 \left(R\left(1 - \sqrt[3]{\frac{\mu}{3}}\right),\; 0\right)\\
L_2 \left(R\left(1 + \sqrt[3]{\frac{\mu}{3}}\right), \; 0\right)\\
L_3 \left(R\left(1 - \frac{5}{12}\mu \right), \; 0 \right)\\
L_4 \left(\frac{R}{2} \cdot \frac{m_1-m_2}{m_1+m_2}, \; \frac{\sqrt{3}}{2}R\right)\\
L_5 \left(\frac{R}{2} \cdot \frac{m_1-m_2}{m_1+m_2}, \; -\frac{\sqrt{3}}{2}R\right)\\
\end{gather}$

где:  

- ㅤ$m_1$ масса большего тела  

- ㅤ$m_2$ масса меньшего тела

- ㅤ$R$ расстояние между телами
"""

# ╔═╡ cfec2969-c486-43e8-aeed-f53157a55563
html"""
<h1 align="center"> Эпизод III<br> <b>Моделирование движения вблизи точек Лагранжа</b></h1><hr>
"""

# ╔═╡ 9a3aedf3-4958-4284-81f5-654f7835fae2
md"""

Движение точки в задаче трёх тел описывается уравнениями движения, описанными во введении. В программе перемещение точки задаётся функцией `move(x, y, px=0, py=0, t=100000, scfale_large)`, где `x` и `y` - обязательные параметры, являющиеся начальными координатами, а `px` и `py` начальные импульсы тела, `t` - общее время движения тела, `scale_large` (если `True`, то показано движение тела относительно системы тел, если `False`, то демонстрируется движение у точки Лагранжа)
"""

# ╔═╡ 7979f857-0489-4435-ad9e-1240c90a3525
begin 
	function euler_iter!(body::Body, delta_t)
		
		push!(body.x, body.x[end] + x_dot(body.x[end], body.y[end], body.px, body.py, body.mu) * delta_t)
		push!(body.y, body.y[end] + y_dot(body.x[end-1], body.y[end], body.px, body.py, body.mu) * delta_t)
		body.px += px_dot(body.x[end], body.y[end], body.px, body.py, body.mu) * delta_t
		body.py += py_dot(body.x[end], body.y[end], body.px, body.py, body.mu) * delta_t
	end
	nothing
end

# ╔═╡ c322b5b7-3a0b-4a91-a833-8ea5fc0a574c
begin
	function leapfrog_iter!(body::Body, x_vel::Float64, y_vel::Float64, px_vel::Float64, py_vel::Float64, delta_t::Float64)
		derivatives = all_dots(body.x[end], body.y[end], body.px, body.py, body.mu, with_impulses=true)
		
		x_vel += derivatives[5] * delta_t
		y_vel += derivatives[6] * delta_t
		px_vel += derivatives[7] * delta_t
		py_vel += derivatives[8] * delta_t
		push!(body.x, body.x[end] + x_vel * delta_t)
		push!(body.y, body.y[end] + y_vel * delta_t)
	
		body.px += px_vel * delta_t
		body.py += py_vel * delta_t
		return (x_vel, y_vel, px_vel, py_vel)
	end
	nothing
end

# ╔═╡ de3a78b5-bcd0-4537-8354-1418f449be4a
# ╠═╡ disabled = true
#=╠═╡
begin
	function move_leapfrog(body::Body, t::Float64=1.0, delta_t::Float64=1e-4)
	    
		x_vel, y_vel, px_vel, py_vel = initial_velocities(body.x[end], body.y[end], body.px, body.py, delta_t)
	    
	    for i in 0:delta_t:t
	        x_vel, y_vel, px_vel, py_vel = leapfrog_iter!(body, x_vel, y_vel, px_vel, py_vel, delta_t)
	    end
	    return (body.x, body.y)
	end
	nothing
end
  ╠═╡ =#

# ╔═╡ 23c2bc64-eb21-4a32-b428-e2d9b325c1c1
function verlet_iter!(body::Body, delta_t::Float64=1e-4)
	derivatives = all_dots(body.x[end], body.y[end], body.px, body.py, body.mu)
	
	push!(body.x, 2 * body.x[end] - body.x[end - 1] + derivatives[5] * delta_t ^ 2)
	push!(body.y, 2 * body.y[end] - body.y[end - 1] + derivatives[6] * delta_t ^ 2)
	
	body.px += derivatives[3] * delta_t
	body.py += derivatives[4] * delta_t
end

# ╔═╡ f09515af-4e9e-48bf-acf2-3acade5c241e
function move_verlet(body::Body, t::Float64, delta_t::Float64=1e-4)
	derivatives = all_dots(body.x[end], body.y[end], body.px, body.py, body.mu)
	push!(body.x, body.x[end] + derivatives[1] * delta_t + derivatives[5] * delta_t^2 / 2)
	push!(body.y, body.y[end] + derivatives[2] * delta_t + derivatives[6] * delta_t^2 / 2)
	for i ∈ 0:delta_t:t
		verlet_iter!(body, delta_t)
	end
	return (body.x, body.y)
end

# ╔═╡ d97e418c-cf1f-4ba2-9887-d5ebfc65576a
function RK2_step!(body::Body, delta_t::Float64=1e-4)
	k1::Vector{Float64} =  collect(all_dots(body.x[end], body.y[end], body.px, body.py, body.mu)[1:end-2])
	k2::Vector{Float64} = collect(initial_velocities(body.x[end], body.y[end], body.px, body.py, body.mu, delta_t))
	x, y, px, py = [body.x[end], body.y[end], body.px, body.py] + (delta_t / 2 * (k1 + k2))
	
	push!(body.x, x)
	push!(body.y, y)
	body.px = px
	body.py = py
end

# ╔═╡ 4d54fb75-8ce1-46a5-958c-79de9572ed49
function move_RK2(body::Body, t::Float64, delta_t::Float64=1e-4)
	for i ∈ 0:delta_t:t
		RK2_step!(body, delta_t)
	end
	return (body.x, body.y)
end

# ╔═╡ 29fce678-6591-4ec6-8d5f-f9b08a8de7ca
function RK4_step!(body::Body, delta_t::Float64=1e-4)
	position::Vector{Float64} = [body.x[end], body.y[end], body.px, body.py]
	k1::Vector{Float64} = collect(all_dots(position..., body.mu, with_second=false))
	k2::Vector{Float64} = collect(all_dots((position .+ delta_t / 2 .* k1)..., body.mu, with_second=false))
	k3::Vector{Float64} = collect(all_dots((position .+ delta_t / 2 .* k2)..., body.mu, with_second=false))
	k4::Vector{Float64} = collect(all_dots((position .+ delta_t .* k3)..., body.mu, with_second=false))
	x, y, px, py = position .+ (delta_t / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4)

	push!(body.x, x)
	push!(body.y, y)
	body.px = px
	body.py = py
	[x, y, px, py]
end

# ╔═╡ 7d5a29c9-501b-468f-9611-d2669ae107da
function move_RK4(body::Body, t::Float64, delta_t::Float64=1e-4)
	for i ∈ 0:delta_t:t
		RK4_step!(body, delta_t)
	end
	return (body.x, body.y)
end

# ╔═╡ 7b828d7b-8505-41c1-a780-1ba84e951872
function stormer_step!(body::Body, delta_t::Float64=1e-4)
	position = [body.x[end], body.y[end], body.px, body.py]
	k1_1 = position + (delta_t / 2) .* p_dots(position..., body.mu)
	k1_2 = k1_1 .+ (q_dots(k1_1..., body.mu) .* delta_t)
	k2 = position + (delta_t / 2) .* (q_dots(k1_1..., body.mu) .+ q_dots(k1_2..., body.mu))
	k3 = k1_1 .+ (delta_t / 2) .* p_dots(k1_2..., body.mu)
	position = [k2[1], k2[2], k3[3], k3[4]]
	change_body!(body, position...)
end

# ╔═╡ 774cf6cb-1a89-4b88-b28e-7c78fc9471a1
function move_stormer(body::Body, t::Float64, delta_t::Float64=1e-4)
for i ∈ 0:delta_t:t
		stormer_step!(body, delta_t)
	end
	return (body.x, body.y)
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
	L1 = move_stormer(Body(x=L_x[1] - 1e-6, y=0.0, px=0.0, py=0.4372524), 1.44)
	Plots.plot(L1[1], L1[2])
	
end

# ╔═╡ 32629dd1-a7fd-45d9-9623-d5665586f489
begin
	L2 = move_stormer(Body(x=L_x[2], y=L_y[2], px=0.0, py=1.27183), 2.15)
	Plots.plot(L2[1], L2[2])
end

# ╔═╡ b8411d90-48e6-403c-ad68-d6671bb53dd6
begin
	L3 = move_stormer(Body(x=L_x[3], y=L_y[3], px=0.0, py=-1.0855), 2.35)
	Plots.plot(L3[1], L3[2])
end

# ╔═╡ 0f50bb6b-55ca-446a-9ee3-18b528b947d4
begin
	L4 = move_stormer(Body(x=L_x[4], y=L_y[4], px=-1.0, py=0.0), 100.0)
	Plots.plot(L4[1], L4[2])
end

# ╔═╡ ec834d95-fba6-4b7f-9ab4-804113dabf1b
begin
	L5 = move_stormer(Body(x=L_x[5], y=L_y[5], px=0.9, py=0.5), 100.0)
	Plots.plot(L5[1], L5[2])
end

# ╔═╡ 93267c29-cc40-4515-a1fb-28c806d4cbb5
begin
	body1 = Body(L_x[1], L_y[1], 0.0, 0.0)
	move_stormer(body1, 30.0)
	Plots.plot(body1.x, body1.y)
	Plots.scatter!([-body1.mu, 1-body1.mu], [0, 0])
end

# ╔═╡ bb78e7c5-53bf-4c3d-8fec-a9f5a6adea5a
begin
	delta_t = 1e-4
	ans_x1, ans_y1 = puancare(body1.x, body1.y, delta_t)
	Plots.scatter(ans_x1, ans_y1)
end

# ╔═╡ bd0514de-cc8d-485a-b0e5-46c8d130988f
begin
	body2 = Body(x = 0.4, y = 0.4, px = 0.334, py = -0.5)
	move_stormer(body2, 100.0)
	
	Plots.plot(body2.x, body2.y)
	Plots.scatter!([-body2.mu, 1-body2.mu], [0, 0])
end

# ╔═╡ d23baf85-ce21-4a93-bc95-a6643cd59ce7
begin
	ans_x2, ans_y2 = puancare(body2.x, body2.y, 1e-4)
	Plots.scatter(ans_x2, ans_y2)
end

# ╔═╡ 8eb16823-767b-4fcf-98da-3da540ca66fd
md"""
## Вывод интеграла Якоби, единственной формы энергии, сохраняющейся в ограниченной задаче трёх тел:

---

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
@fastmath function zero_velocity(body::Body; t::Float64=1000.0, ϵ::Float64=0.01)
	ans_xx::Array{Float64}, ans_yy::Array{Float64} = [], []
	ϵ::Float64 = 0.01
	const_h = h_const(body_cords(body)...)
	for i in -2.0:0.01:2.0
	    for j in -2.0:0.01:2.0
	        if abs(-2u(i, j) + const_h) <= ϵ
	            append!(ans_xx, i)
	            append!(ans_yy, j)
	        end
	    end
	end
	a, b = move_stormer(body, 1000.0)
	scatter(ans_xx, ans_yy)
	plot!(a, b)
end

# ╔═╡ 919be67a-1988-4cf7-a660-6db6183fef1c
zero_velocity(Body(x=0.2, y=0.4, px=0.0, py=0.0))

# ╔═╡ 8a948b7e-8520-45b4-93ef-3add2d0c50da
zero_velocity(Body(x=0.5, y=0.5, px=-0.2, py=0.2))

# ╔═╡ 5a001069-c770-4450-b09f-f0b4f1f7f5c1
# ╠═╡ disabled = true
#=╠═╡
zero_velocity(Body(x=L_x[5], y=L_y[5], px=0.9, py=0.5))
  ╠═╡ =#

# ╔═╡ 5bec7dcf-5d47-4674-9cf9-7dabb4c671bf
md"""
Воспользуемся интегралом Якоби для определения точности численного интегрирования
"""

# ╔═╡ 1f65f3e5-432e-4c42-b194-6b16f429ce38
# ╠═╡ disabled = true
#=╠═╡
function accuracy_leapfrog(body::Body; t::Float64=10.0, delta_t::Float64=1e-4)
	h_const0 = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
	ans_h = []
	x_vel, y_vel, px_vel, py_vel = initial_velocities(body.x[end], body.y[end], body.px, body.py, delta_t)
	for i ∈ 1:delta_t:t
		x_vel, y_vel = leapfrog_iter!(body, x_vel, y_vel, px_vel, py_vel, delta_t)
		h_new = h_const(body.x[end], body.y[end], body.px, body.py, body.mu) 
		push!(ans_h, (abs(h_new) - abs(h_const0)) / abs(h_const0) * 100)
	end
	return (ans_h, 1:delta_t:t)
	
end
  ╠═╡ =#

# ╔═╡ 024950bd-c4e0-454d-bf51-4c2bac424f90
#=╠═╡
begin
	ak1, tk1 = accuracy_leapfrog(Body(L_x[1], L_y[1], 0.0, 0.0), t=200.0)
	scatter(tk1, ak1, markersize=0.000001)
end
  ╠═╡ =#

# ╔═╡ 86e576f7-32fa-428c-b921-1d4652c009aa
function accuracy_euler(body::Body; t::Float64=10.0, delta_t::Float64=1e-4)
	h_const0 = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
	ans_h = []
	for i ∈ 1:delta_t:t
		euler_iter!(body, delta_t)
		h_new = h_const(body.x[end], body.y[end], body.px, body.py, body.mu) 
		push!(ans_h, (abs(h_new) - abs(h_const0)) / abs(h_const0) * 100)
	end
	return (ans_h, 1:delta_t:t)
	
end

# ╔═╡ 68f05f3f-01e7-49af-af05-fa2a8f09fe4d
# ╠═╡ disabled = true
#=╠═╡
function accuracy_verlet(body::Body, t::Float64, delta_t::Float64=1e-4)
	derivatives = all_dots(body.x[end], body.y[end], body.px, body.py, body.mu)
	h_const0 = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
	ans = [h_const0]
	
	push!(body.x, body.x[end] + derivatives[1] * delta_t + derivatives[5] * delta_t^2 / 2)
	push!(body.y, body.y[end] + derivatives[2] * delta_t + derivatives[6] * delta_t^2 / 2)
	push!(ans, h_const(body.x[end], body.y[end], body.px, body.py, body.mu))
	for i ∈ 0:delta_t:t
		verlet_iter!(body, delta_t)
		h_new = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
		push!(ans, (h_new - h_const0) / h_const0 * 100) 

	end
	return (ans, 0:delta_t:t)
end
  ╠═╡ =#

# ╔═╡ 1127fe48-8f43-40be-a861-00cc641ec3fb
# ╠═╡ disabled = true
#=╠═╡
begin
	ak3, tk3 = accuracy_verlet(Body(L_x[1], L_y[1], 0.0, 0.0), 200.0)
	scatter(tk3, ak3, markersize=0.00001)
end
  ╠═╡ =#

# ╔═╡ 08f3874a-8a2a-4c16-a9be-749beee1547b
# ╠═╡ disabled = true
#=╠═╡
begin
	ak2, tk2 = accuracy_euler(Body(L_x[1], 0.4, 0.0, 0.0), t=200.0)
	scatter(tk2, ak2, markersize=0.00001)
end
  ╠═╡ =#

# ╔═╡ ab40a600-e4ca-4927-a8cf-85d98cc9e7ca
# ╠═╡ disabled = true
#=╠═╡
function accuracy_RK2(body::Body; t::Float64, delta_t::Float64=1e-4)
	h_const0 = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
	ans = []
	for i ∈ 0:delta_t:t
		RK4_step!(body, delta_t)
		h_new = h_const(body.x[end], body.y[end], body.px, body.py, body.mu)
		push!(ans, (h_new - h_const0) / h_const0 * 100) 

	end
	return (ans, 0:delta_t:t)
end
  ╠═╡ =#

# ╔═╡ defebd59-ba85-4b8b-a4e4-00e364b87f6e
#=╠═╡
begin
	ak2, tk2 = accuracy_RK2(Body(L_x[1], 0.4, 0.4, 0.0), t=100.0)
	scatter(tk2, ak2, markersize=0.00001)
end
  ╠═╡ =#

# ╔═╡ 30358f2c-0329-4b59-82e4-f208e8596f25
md"""
Тест при большом расстоянии третьего тела от объектов 
"""

# ╔═╡ c361a162-30ee-47b8-88f4-2e59225da70c
begin
	ɼ1, ɼ2 = move_stormer(Body(x=1000.0, y=10000000.0, px=1000000.0, py=10000000.0), 200.0)
	plot(ɼ1, ɼ2)
end

# ╔═╡ 56d02fea-a060-4634-abe9-082c706fc0f3
md"""
Построение анимации в сидерической системе отсчёта
"""

# ╔═╡ 1c0ad33c-2490-4922-9c61-8b5595cef356
@fastmath function animate_non_rotating(body::Body; t::Float64=10.0, delta_t::Float64=1e-4, frames::Int64 = 200, xlim::Float64=2.0, ylim::Float64=2.0, name::String="main.gif", fps::Int64=60)
	ro1, ro2 = -body.mu, 1 - body.mu
	psi = 0
	circle_dots = range(0, 2π, length = 100)

		
	delta_t = 1e-4
    x_vel, y_vel = initial_velocities(body.x[end], body.y[end], body.px, body.py, delta_t)
	points = [-body.mu, 1-body.mu], [0, 0]
	scatter(points)
	
	anim = @animate for i in 0:1:(t / delta_t)
        x_vel, y_vel = leapfrog_iter!(body, x_vel, y_vel, delta_t)
		psi += delta_t 
		
		
		if i % frames == 0
			x1, y1 = ro1 * cos(psi), ro1 * sin(psi)
		x2, y2 = ro2 * cos(psi), ro2 * sin(psi)
			x_1 = body.x[end] * cos(psi) - body.y[end] * sin(psi)
			y_1 = body.x[end] * sin(psi) + body.y[end] * cos(psi)
			scatter([x1, x2], [y1, y2], xlim=(-xlim, xlim), ylim=(-ylim, ylim), color="black")
			scatter!([x_1], [y_1], color="magenta")
			plot!(@. sin(circle_dots) * ro1, @. cos(circle_dots) * ro1) 
			plot!(@. sin(circle_dots) * ro2, @. cos(circle_dots) * ro2) 
		end
	end every frames
	gif(anim, name, fps=fps)
end

# ╔═╡ 984d2d1b-8333-4e88-ba23-d5fdcc7927a4
md"""
**Проверка теоремы Лиувилля о сохранение фазового объёма**
"""

# ╔═╡ efd1bb93-bf79-4aef-aa2c-a3622750d632
function move_some(bodies::Array{Body}, t:: Float64=5.0; delta_t::Float64=1e-4, move::Function=move_RK4)
	for i ∈ bodies
		move(i, t, delta_t)
	end
end

# ╔═╡ 292ffba1-ef46-4e66-88ce-783618467592
function polygon_area(dots::Array)
	"""Функция нахождения площади выпуклой оболочки"""
	x0 = dots[1, 1]
	y0 = dots[1, 2]
	sq = 0
	
	for i ∈ 2:size(dots)[1]
		sq += (dots[i, 1] - dots[i - 1, 1]) * (dots[i, 2] + dots[i - 1, 2])
	end
	sq += (dots[end, 1] - x0) * (dots[end, 2] + y0)
	return abs(sq / 2)
end

# ╔═╡ ce716daa-c750-41a9-86c3-6b63f88903c6
polygon_area([1.0 1.0; 4.0 1.0; 2.0 5.0])

# ╔═╡ 10b0371e-fc41-4696-bbc7-440b26236274
function hull(bodies::Array{Body})

end

# ╔═╡ f311fb12-9555-4718-b35d-d5be82352cdc
function liouville_x(x, y, px, py, mu=0.2; dx=1e-3, dpx=1e-3, t=10.0, step_t=1.0)
	body_lst = [Body(x, y, px, py, mu), Body(x+dx, y, px, py,mu), Body(x, y, px+dpx, py, mu)]
	ans = [dpx * dx]
	for i ∈ 0:step_t:t
		for body ∈ body_lst
			move_RK2(body, step_t)
		end
		push!(ans, (body_lst[2].x[end] - body_lst[1].x[end]) * (body_lst[3].px - body_lst[1].px))
		
	end
	println(ans)
	println(collect(0:step_t:t))
	return (collect(0:step_t:t), ans)
end

# ╔═╡ 45ebe637-d840-4f7b-8e64-04f8d3fed5ef
md"""
_Интеграл Якоби_

$C_j = 2U(x, y) - (\dot x^2 + \dot y^2)$

Произведём замену $\dot x = p_x + y$, $\dot y = p_y - x$
Посдставим значения в интеграл

$C_j = 2U(x, y) - ((p_x + y)^2 + (p_y - x)^2)$

Пусть 

$\dot x = p_x + y = 0$ 
"""

# ╔═╡ 4f4728fe-e264-4231-88e0-a89258d79b61
function Jacobi_move(x::Float64, y::Float64, Cⱼ::Float64, μ::Float64=0.2)::Body
	px = - y
	py = sqrt(abs(2*u(x, y, μ) + Cⱼ)) + x
	return Body(x, y, px, py, μ)
end

# ╔═╡ 8f9b8829-ddd9-448d-b728-374d995ab59d
# ╠═╡ disabled = true
#=╠═╡
begin
	const LAGRANGE_1_2 = Lagrange(1/82.27)
	bodyₙ = Jacobi_move(1-1/82, 0.0, -0.8668, 1/82.27)
	plot(move_RK2(bodyₙ, 10.0)...)
	scatter!([-0.5, 0.5], [0, 0])
end
  ╠═╡ =#

# ╔═╡ d5844f01-7cd0-4105-b72e-1dec9695cc54
md"""
# **Поиск периодических орбит**
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
ProgressBars = "49802e3a-d2f1-5c88-81d8-b72133a6f568"

[compat]
Parameters = "~0.12.3"
Plots = "~1.39.0"
ProgressBars = "~1.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "5628a1d30e731d4c623d3ef79c5beed952aefdf3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d73afa4a2bb9de56077242d98cf763074ab9a970"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1596bab77f4f073a14c62424283e7ebff3072eca"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.9+1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "cb56ccdd481c0dd7f975ad2b3b62d9eda088f7e2"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.14"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bbb5c2115d63c2f1451cb70e5ef75e8fe4707019"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.22+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressBars]]
deps = ["Printf"]
git-tree-sha1 = "b437cdb0385ed38312d91d9c00c20f3798b30256"
uuid = "49802e3a-d2f1-5c88-81d8-b72133a6f568"
version = "1.5.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "364898e8f13f7eaaceec55fd3d08680498c0aa6e"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.4.2+3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═f654666c-36bb-4854-ae63-3f0a1bb0343f
# ╟─d3f35948-a39c-462c-9412-8bbefa334fe2
# ╟─a3e811d7-0a06-44a6-a9db-0d21983b212a
# ╠═f997a055-a870-430a-8aca-5060c2fd16cb
# ╠═d3f8fd38-ba29-4234-a316-724cfb5ff183
# ╠═56257e88-88f1-4c52-a9b4-382f6de339f6
# ╟─c1590f21-e4eb-48cb-8c0d-fb98fbfed059
# ╠═fc385955-4924-4326-ab8e-d37ecd5f65e4
# ╠═41788216-b144-4460-bded-e87992c3d57a
# ╠═e41e33a5-212f-403f-b51c-ae94dfab293b
# ╟─9f7c73f7-226e-4c04-9741-834f46531deb
# ╟─8ff68780-6554-4873-9abb-aa118e116f62
# ╠═2a4aa267-796c-4090-ad20-3ce806868726
# ╠═7befb0a2-e81f-4d8b-bd09-4f253dc36325
# ╟─1e8db92f-d780-4052-9381-87d0719659cb
# ╟─ad124e1c-e3bf-47f7-b7a5-f1d97e00a14c
# ╟─cfec2969-c486-43e8-aeed-f53157a55563
# ╟─9a3aedf3-4958-4284-81f5-654f7835fae2
# ╠═7979f857-0489-4435-ad9e-1240c90a3525
# ╠═c322b5b7-3a0b-4a91-a833-8ea5fc0a574c
# ╠═de3a78b5-bcd0-4537-8354-1418f449be4a
# ╠═23c2bc64-eb21-4a32-b428-e2d9b325c1c1
# ╠═f09515af-4e9e-48bf-acf2-3acade5c241e
# ╠═d97e418c-cf1f-4ba2-9887-d5ebfc65576a
# ╠═4d54fb75-8ce1-46a5-958c-79de9572ed49
# ╠═29fce678-6591-4ec6-8d5f-f9b08a8de7ca
# ╠═7d5a29c9-501b-468f-9611-d2669ae107da
# ╠═7b828d7b-8505-41c1-a780-1ba84e951872
# ╠═774cf6cb-1a89-4b88-b28e-7c78fc9471a1
# ╟─911f528b-816a-4c6e-9538-2530309f70c6
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
# ╠═919be67a-1988-4cf7-a660-6db6183fef1c
# ╠═8a948b7e-8520-45b4-93ef-3add2d0c50da
# ╠═5a001069-c770-4450-b09f-f0b4f1f7f5c1
# ╟─5bec7dcf-5d47-4674-9cf9-7dabb4c671bf
# ╠═1f65f3e5-432e-4c42-b194-6b16f429ce38
# ╠═024950bd-c4e0-454d-bf51-4c2bac424f90
# ╠═86e576f7-32fa-428c-b921-1d4652c009aa
# ╠═68f05f3f-01e7-49af-af05-fa2a8f09fe4d
# ╠═1127fe48-8f43-40be-a861-00cc641ec3fb
# ╠═08f3874a-8a2a-4c16-a9be-749beee1547b
# ╠═ab40a600-e4ca-4927-a8cf-85d98cc9e7ca
# ╠═defebd59-ba85-4b8b-a4e4-00e364b87f6e
# ╟─30358f2c-0329-4b59-82e4-f208e8596f25
# ╠═c361a162-30ee-47b8-88f4-2e59225da70c
# ╟─56d02fea-a060-4634-abe9-082c706fc0f3
# ╟─1c0ad33c-2490-4922-9c61-8b5595cef356
# ╟─984d2d1b-8333-4e88-ba23-d5fdcc7927a4
# ╟─efd1bb93-bf79-4aef-aa2c-a3622750d632
# ╟─292ffba1-ef46-4e66-88ce-783618467592
# ╟─ce716daa-c750-41a9-86c3-6b63f88903c6
# ╟─10b0371e-fc41-4696-bbc7-440b26236274
# ╟─f311fb12-9555-4718-b35d-d5be82352cdc
# ╟─45ebe637-d840-4f7b-8e64-04f8d3fed5ef
# ╠═4f4728fe-e264-4231-88e0-a89258d79b61
# ╠═8f9b8829-ddd9-448d-b728-374d995ab59d
# ╟─d5844f01-7cd0-4105-b72e-1dec9695cc54
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
