from django.urls import reverse, resolve
from django.shortcuts import redirect


class RedirectToLoginMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        url_name = resolve(request.path_info).url_name

        if (
                request.path.startswith('/api/')
                and not request.user.is_authenticated
                and 'Authorization' not in request.headers
                and url_name != 'login'
                and url_name != 'token_obtain_pair'
                and url_name != 'token_refresh'
        ):
            return redirect(reverse('api_login'))

        response = self.get_response(request)
        return response