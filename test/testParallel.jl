@testset "ParallelBackend" begin
    @testset "detect_resources returns NamedTuple" begin
        res = detect_resources()
        @test res isa NamedTuple
        @test haskey(res, :nthreads)
        @test haskey(res, :nworkers)
        @test haskey(res, :pids)
        @test haskey(res, :has_cuda)
        @test haskey(res, :gpu_count)
        @test haskey(res, :has_darrays)
        @test haskey(res, :has_parallelops)
        # Invariants
        @test res.nthreads >= 1
        @test res.nworkers >= 0
        @test res.gpu_count >= 0
        @test res.has_darrays isa Bool
        @test res.has_parallelops isa Bool
    end

    @testset "select_backend defaults to serial/threads" begin
        b = select_backend()
        @test b isa ParallelBackend
        @test b.kind in (:serial, :threads, :distributed, :gpu)
        @test b.gpu isa Bool
        @test b.nthreads >= 1
        @test b.nworkers >= 0
        @test b.fft_threads >= 1
        # In a CI box without workers, kind should not be :distributed.
        if Distributed.nworkers() < 2
            @test b.kind != :distributed
        end
    end

    @testset "select_backend: explicit kind=serial" begin
        b = select_backend(kind=:serial)
        @test b.kind === :serial
        @test is_serial(b)
        @test !is_threads(b)
        @test !is_distributed(b)
        @test !is_gpu(b)
        @test !is_parallel(b)
    end

    @testset "select_backend: explicit kind=threads (or falls back)" begin
        # When only 1 thread is available, this falls back to serial with a warning.
        b = select_backend(kind=:threads)
        @test b.kind in (:threads, :serial)
    end

    @testset "select_backend: unknown kind → ArgumentError" begin
        @test_throws ArgumentError select_backend(kind=:bogus)
    end

    @testset "convenience predicates" begin
        s = select_backend(kind=:serial)
        @test is_serial(s)
        @test !is_parallel(s)
        @test_throws MethodError is_serial(:not_a_backend)  # predicate is only for ParallelBackend
    end

    @testset "device_array / local_array are functions" begin
        b = select_backend(kind=:serial)
        @test b.device_array isa Function
        @test b.local_array isa Function
        @test b.release! isa Function
        @test b.distribute_array isa Function
    end

    @testset "serial backend: device_array / local_array are identity" begin
        b = select_backend(kind=:serial)
        A = [1.0, 2.0, 3.0]
        @test b.device_array(A) === A
        @test b.local_array(A) === A
        @test b.distribute_array(A) === A
        # release! must return nothing
        @test b.release!(A) === nothing
    end

    @testset "serial backend: fft_threads == 1" begin
        b = select_backend(kind=:serial)
        @test b.fft_threads == 1
    end

    @testset "set_backend_fft_threads! pins FFTW thread count" begin
        b = select_backend(kind=:serial)
        pinned = set_backend_fft_threads!(b)
        @test pinned == 1
        # After pinning, FFTW should report the same thread count
        @test FFTW.get_num_threads() == 1
    end

    @testset "Base.show for ParallelBackend" begin
        b = select_backend(kind=:serial)
        s = sprint(show, b)
        # `kind` is rendered as the first positional arg, not as a `kind=` kw.
        # The remaining fields use `name=value` form (see backend.jl show method).
        @test occursin("ParallelBackend", s)
        @test occursin("serial", s)
        @test occursin("gpu=", s)
        @test occursin("nthreads=", s)
        @test occursin("nworkers=", s)
    end
end

@testset "Parallel operations" begin
    b = select_backend(kind=:serial)

    @testset "to_device / to_host identity on serial backend" begin
        A = rand(3, 3)
        @test to_device(A, b) === A
        @test to_host(A, b) === A
    end

    @testset "distribute identity on serial backend" begin
        A = rand(5)
        @test distribute(A, b) === A
    end

    @testset "release! is no-op on serial backend" begin
        A = rand(2, 2)
        @test release!(A, b) === nothing
        # The original array must still be usable
        @test sum(A) isa Real
    end

    @testset "parallel_sum / maximum / minimum forward to Base" begin
        A = [1.0, 2.0, 3.0, 4.0]
        @test parallel_sum(A, b) == sum(A)
        @test parallel_maximum(A, b) == maximum(A)
        @test parallel_minimum(A, b) == minimum(A)
    end

    @testset "parallel_findmax returns (value, index)" begin
        A = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0]
        v, idx = parallel_findmax(A, b)
        @test v == 9.0
        @test idx == 6
    end

    @testset "parallel_quantile returns expected quantile" begin
        values  = [1.0, 2.0, 3.0, 4.0, 5.0]
        weights = [1.0, 1.0, 1.0, 1.0, 1.0]
        q = parallel_quantile(values, weights, [0.5], b)
        # Linear-interpolated weighted quantile at 0.5 of uniform weights
        @test length(q) == 1
        @test q[1] ≈ 3.0 atol=1e-9
    end

    @testset "parallel_broadcast! is in-place" begin
        dst = zeros(4)
        src = [1.0, 2.0, 3.0, 4.0]
        parallel_broadcast!(identity, dst, (src,), b)
        @test dst == src
    end

    @testset "parallel_sumprod = dot product" begin
        A = [1.0, 2.0, 3.0]
        B = [4.0, 5.0, 6.0]
        @test parallel_sumprod(A, B, b) ≈ sum(A .* B)
    end
end

@testset "Parallel backend kind invariants" begin
    @testset "kind field accepts only known symbols" begin
        b = select_backend(kind=:serial)
        @test b.kind === :serial
        b2 = select_backend(kind=:auto)
        @test b2.kind in (:serial, :threads, :distributed, :gpu)
    end

    @testset "gpu flag tracks kind" begin
        for kind in (:serial, :threads)
            b = select_backend(kind=kind)
            @test b.gpu === false
            @test !is_gpu(b)
        end
    end

    @testset "pids is empty for non-distributed backends" begin
        b = select_backend(kind=:serial)
        @test b.pids == Int[]
        b2 = select_backend(kind=:threads)
        @test b2.pids == Int[]
    end
end
