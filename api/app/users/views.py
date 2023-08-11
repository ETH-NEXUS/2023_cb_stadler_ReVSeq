import json
from django.contrib.auth import authenticate, login, logout
from django.http import JsonResponse
from django.utils.translation import gettext as _
from django.core.handlers.wsgi import WSGIRequest
from django.views.generic import View
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import ensure_csrf_cookie
from rest_framework.response import Response
from django.contrib.auth import get_user_model
from .serializers import UserSerializer

from rest_framework import viewsets
import logging

logger = logging.getLogger(__name__)


class UserViewSet(viewsets.ModelViewSet):
    queryset = get_user_model().objects.all()
    serializer_class = UserSerializer
    lookup_field = "pk"

    def retrieve(self, request, *args, **kwargs):
        instance = self.request.user
        serializer = self.get_serializer(instance)
        return Response(serializer.data)


class CsrfCookieView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request: WSGIRequest, *args, **kwargs):
        return JsonResponse({"details": _("CSRF cookie set")})


class LoginView(View):
    def post(self, request: WSGIRequest, *args, **kwargs):
        logger.debug("Start login")
        data = json.loads(request.body)
        username = data.get("username")
        password = data.get("password")
        print(username, password)

        if username is None or password is None:
            return JsonResponse(
                {"detail": _("Please provide username and password.")}, status=400
            )
        print("Before authenticate")
        user = authenticate(username=username, password=password)

        if user is None:
            return JsonResponse({"detail": _("Invalid credentials.")}, status=400)
        print("After authenticate")
        login(request, user)
        print("After login")
        logger.debug("End login")
        return JsonResponse({"detail": _("Successfully logged in.")})


class LogoutView(View):
    def get(self, request: WSGIRequest, *args, **kwargs):
        if not request.user.is_authenticated:
            return JsonResponse({"detail": _("You're not logged in.")}, status=400)

        logout(request)
        return JsonResponse({"detail": _("Successfully logged out.")})
